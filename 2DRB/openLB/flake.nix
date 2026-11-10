{
  description = "OpenLB";

  inputs = {
    nixpkgs.url = github:NixOS/nixpkgs/nixos-25.11;
  };

  outputs = { self, nixpkgs, ... }: let
    system = "x86_64-linux";
    lib = nixpkgs.lib;

    cudaMpiOverlay = (final: prev: {
      ucx-cuda = (prev.ucx.override {
        enableCuda = true;
      }).overrideAttrs (old: {
        # UCX must be built multi-thread-capable, else OpenMPI cannot select the UCX PML
        # under OpenLB's required MPI_THREAD_MULTIPLE and silently falls back to OB1+smcuda
        # -> intra-node GPU halos host-stage over PCIe (NVLink cuda_ipc never used).
        configureFlags = (old.configureFlags or []) ++ [ "--enable-mt" ];
      });

      openmpi-cuda = (prev.openmpi.override {
        cudaSupport = true;
        ucx = final.ucx-cuda;
      }).overrideAttrs (old:
        let
          cuda-merged = final.symlinkJoin {
            name = "cuda-merged";
            paths = with final.cudaPackages; [
              cuda_nvcc
              cuda_cudart
              cuda_cccl
            ];
          };
        in {
          nativeBuildInputs = (old.nativeBuildInputs or []) ++ [ final.cudaPackages.cuda_nvcc ];
          buildInputs = (old.buildInputs or []) ++ [ cuda-merged ];
          configureFlags = (old.configureFlags or []) ++ [
            "--with-cuda=${cuda-merged}"
            "--with-cuda-libdir=${cuda-merged}/lib/stubs"
          ];
        });

      fixedPythonPackages = prev.python3Packages.override {
        overrides = pyFinal: pyPrev: {
          mpi4py = pyPrev.mpi4py.overrideAttrs (old: {
            doCheck = false;
            checkPhase = "true";
            installCheckPhase = "true";
            doInstallCheck = false;
          });
        };
      };

      adios2-cuda = (prev.adios2.override {
        mpi = final.openmpi-cuda;
        mpiSupport = true;
        ucx = final.ucx-cuda;
        python3Packages = final.fixedPythonPackages;
        pythonSupport = true;
      });
    });

    pkgs = import nixpkgs {
      inherit system;
      config.allowUnfree = true;
      overlays = [ cudaMpiOverlay ];
    };

    # Tools present in every dev shell.
    baseTools = with pkgs; [
      gnumake
      git
      gdb
      valgrind
      clang-tools
      gnuplot

      python3
      uv

    # meshing
      kahip
      metis
    ];

    cudaTools = with pkgs; [
      cudaPackages.cudatoolkit
      cudaPackages.cuda_cudart
    ];

    # Link-time CUDA driver stub so -lcuda resolves on driverless runners
    cudaStubLib = "${pkgs.cudaPackages.cuda_cudart}/lib/stubs";

    # Target microarchitecture, resolved on the machine the shell runs on: default
    # to the local CPU but let an HPC build whose login node differs from the
    # compute node pin it via OLB_CXX_MARCH=<arch>. Compiles to bash ${OLB_CXX_MARCH:-...}.
    marchExpr = arch: "\${OLB_CXX_MARCH:-${arch}}";

    # Compute capability of the first local GPU (override OLB_CXX_MARCH's GPU analogue
    # CUDA_ARCH when in doubt). Detected at shell entry, like the CPU arch above.
    cudaArchLine = ''export CUDA_ARCH=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -n 1 | grep -o [0-9] | tr -d '\n')'';

    # HIP-on-CUDA (HIP_PLATFORM=nvidia) needs the nvidia_detail HIP->CUDA headers,
    # which nixpkgs' (AMD-only) clr does not ship. They live in the separate
    # ROCm/hipother repo at hipnv/include/hip/nvidia_detail; pin it to the clr
    # version so the API surface matches.
    hipother = pkgs.fetchFromGitHub {
      owner = "ROCm";
      repo = "hipother";
      rev = "rocm-6.4.3";
      hash = "sha256-jxMiG22fNANO9GBSWlvXH7abZHs3kyQffwVP6DKhPtQ=";
    };

    envDefaults = {
      compiler = "gcc";      # "gcc" | "clang" | "nvcc" (nvcc = sole CUDA driver)
      parallel = "none";     # "none" | "mpi" | "omp" | "hybrid"
      cuda = false;
      adios = false;
      nlopt = false;
      vtk = false;
      profiling = false;
      simd = false;
      arch = "native";
    };

    # Concise, feature-descriptive shell name, e.g. gcc-mpi-cuda-adios. nvcc (the
    # sole-driver CUDA compiler) is simply named "cuda".
    envName = cfg: let
      c = envDefaults // cfg;
      isNvcc = c.compiler == "nvcc";
    in lib.concatStringsSep "-" (
      [ (if isNvcc then "cuda" else c.compiler) ]
      ++ lib.optional (c.parallel == "mpi") "mpi"
      ++ lib.optional (c.parallel == "omp") "omp"
      ++ lib.optional (c.parallel == "hybrid") "hybrid"
      ++ lib.optional c.simd "simd"
      ++ lib.optional (c.cuda && !isNvcc) "cuda"
      ++ lib.optional c.nlopt "nlopt"
      ++ lib.optional c.adios "adios"
      ++ lib.optional c.vtk "vtk"
      ++ lib.optional c.profiling "prof"
    );

    # Compose a dev shell from a feature selection; every axis is derived here so
    # the shell declarations stay one line of intent.
    mkEnv = cfg: let
      c = envDefaults // cfg;

      isClang = c.compiler == "clang";
      isNvcc = c.compiler == "nvcc";   # nvcc as the sole compiler driver
      cudaOn = c.cuda || isNvcc;

      mpi = c.parallel == "mpi" || c.parallel == "hybrid";
      omp = c.parallel == "omp" || c.parallel == "hybrid";

      cc = if isNvcc then "nvcc" else if isClang then "clang" else "gcc";
      cxx = if isNvcc then "nvcc" else if isClang then "clang++" else if mpi then "mpic++" else "g++";
      stdenv = if isClang then pkgs.llvmPackages_19.libcxxStdenv else pkgs.stdenvNoCC;

      platforms = [ "CPU_SISD" ]
        ++ lib.optional c.simd "CPU_SIMD"
        ++ lib.optional cudaOn "GPU_CUDA";
      parallelMode = { none = "NONE"; mpi = "MPI"; omp = "OMP"; hybrid = "HYBRID"; }.${c.parallel};

      features = lib.optional c.adios "ADIOS2"
        ++ lib.optional c.nlopt "NLOPT"
        ++ lib.optional c.vtk "VTK"
        ++ lib.optional c.profiling "PROFILING";

      buildInputs = baseTools
        ++ lib.optional (!isClang) (if cudaOn then pkgs.cudaPackages.backendStdenv.cc else pkgs.gcc13)   # CUDA: nvcc's own host gcc, else gcc13
        ++ lib.optionals isClang (with pkgs; [ clang_19 llvmPackages_19.bintools-unwrapped llvmPackages_19.openmp ])
        ++ lib.optional mpi (if cudaOn then pkgs.openmpi-cuda else pkgs.openmpi)
        ++ lib.optionals cudaOn cudaTools
        ++ lib.optional c.adios pkgs.adios2-cuda
        ++ lib.optional c.nlopt pkgs.nlopt
        ++ lib.optional c.vtk pkgs.vtk;

      useRpath = (cudaOn && mpi) || c.adios || c.vtk;
      rpathDirs = lib.optionals (cudaOn && mpi) [ "${pkgs.openmpi-cuda}/lib" "${pkgs.cudaPackages.cuda_cudart}/lib" ]
        ++ lib.optional c.adios "${pkgs.adios2-cuda}/lib"
        ++ lib.optional c.vtk "${pkgs.vtk}/lib"
        ++ [ "." "/run/opengl-driver/lib" ];

      # -O level honours OLB_CXX_OPT_LEVEL (default 3); CI compile-check jobs set 0 for speed
      cxxflags = "-O\${OLB_CXX_OPT_LEVEL:-3} -Wall -march=$ARCH -mtune=$ARCH -std=c++20"
        + lib.optionalString isClang " -Wno-vla-cxx-extension"
        + lib.optionalString c.adios " $(adios2-config --mpi --cxx-flags)"
        + lib.optionalString c.nlopt " -I${pkgs.nlopt}/include/nlopt"
        + lib.optionalString c.vtk " -I${pkgs.vtk}/include/vtk";

      ldExtra = lib.optionalString useRpath "$RPATH"
        + lib.optionalString c.adios " $(adios2-config --mpi --cxx-libs)"
        + lib.optionalString c.nlopt " -L${pkgs.nlopt}/lib -lnlopt"
        + lib.optionalString c.vtk " -L${pkgs.vtk}/lib";

      ldLibPath = lib.concatStringsSep ":" (lib.optional mpi "." ++ [ "/run/opengl-driver/lib" ]);

      # nvcc is itself the CUDA compiler, so it needs no host CUDA_CXX split.
      cudaBlock = lib.concatStringsSep "\n" (
        lib.optionals (!isNvcc) [ "export CUDA_CXX=nvcc" ''export CUDA_CXXFLAGS="-O3 -std=c++20${lib.optionalString c.adios " $(adios2-config --mpi --cxx-flags)"}${lib.optionalString c.vtk " -I${pkgs.vtk}/include/vtk"}${lib.optionalString omp " --forward-unknown-to-host-compiler"}"'' ]
        ++ [ ''export CUDA_LDFLAGS="-L/run/opengl-driver/lib -L${cudaStubLib}"'' cudaArchLine ''export FLOATING_POINT_TYPE=float'' ]
      );

      shellHook = lib.concatStringsSep "\n" (
        # nvcc rejects the host -march flags, so it gets no CXXFLAGS line (the
        # Makefile supplies its defaults), like the original CUDA environment.
        lib.optional (!isNvcc) ''export ARCH="${marchExpr c.arch}"''
        ++ [
          "export CXX=${cxx}"
          "export CC=${cc}"
        ]
        ++ lib.optional (!isNvcc) ''export CXXFLAGS="${cxxflags}"''
        ++ [ "export PARALLEL_MODE=${parallelMode}" ]
        # Force the UCX PML so intra-node GPU halos use NVLink cuda_ipc. Gated on cudaOn:
        # only the CUDA shells get the --enable-mt ucx-cuda that MPI_THREAD_MULTIPLE needs;
        # forcing it on non-CUDA MPI's stock single-threaded UCX aborts with no fallback.
        ++ lib.optional (mpi && cudaOn) ''export OMPI_MCA_pml=ucx''
        ++ lib.optional omp ''export OMPFLAGS="-fopenmp"''
        ++ [ ''export PLATFORMS="${lib.concatStringsSep " " platforms}"'' ]
        # uv-provisioned standalone CPython (meshing tests, profiling tooling) has no
        # RPATH to the nix C++ stdlib, so its numpy/numba/pyvista wheels fail to load
        # libstdc++.so.6. Put it on the loader path for every shell.
        ++ [ ''export LD_LIBRARY_PATH="''${LD_LIBRARY_PATH:+$LD_LIBRARY_PATH:}${pkgs.stdenv.cc.cc.lib}/lib"'' ]
        ++ lib.optional cudaOn cudaBlock
        ++ lib.optional useRpath ''export RPATH="-Wl,-rpath,${lib.concatStringsSep ":" rpathDirs}"''
        ++ lib.optional (features != []) ''export FEATURES="${lib.concatStringsSep " " features}"''
        ++ lib.optional (ldExtra != "") ''export LDFLAGS="$LDFLAGS ${ldExtra}"''
        ++ lib.optional cudaOn ''export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:${ldLibPath}"''
      );
    in (pkgs.mkShell.override { inherit stdenv; }) {
      name = "openlb-${envName cfg}";
      IN_OLB_NIX_SHELL = "1";
      inherit buildInputs shellHook;
    };

    # Curated mainstream feature combinations; the builder derives the rest.
    bases = [
      { }                                                      # gcc
      { simd = true; }                                         # gcc-simd
      { parallel = "mpi"; }                                    # gcc-mpi
      { parallel = "omp"; }                                    # gcc-omp
      { compiler = "clang"; }                                  # clang
      { compiler = "nvcc"; }                                   # cuda (nvcc sole driver)
      { cuda = true; }                                         # gcc-cuda
      { parallel = "mpi"; cuda = true; }                       # gcc-mpi-cuda
      { parallel = "mpi"; cuda = true; adios = true; }         # gcc-mpi-cuda-adios
      { parallel = "mpi"; cuda = true; vtk = true; }           # gcc-mpi-cuda-vtk
      { parallel = "hybrid"; cuda = true; vtk = true; }        # gcc-hybrid-cuda-vtk
      { parallel = "mpi"; nlopt = true; }                      # gcc-mpi-nlopt
    ];

    # Every base also gets a profiling permutation.
    withProf = base: [ base ] ++ lib.optional (!(base.profiling or false)) (base // { profiling = true; });

    generated = lib.listToAttrs (
      map (cfg: lib.nameValuePair (envName cfg) (mkEnv cfg)) (lib.concatMap withProf bases)
    );

    # FEATURES export line from a list, adding PROFILING for prof permutations.
    featuresLine = feats: profiling:
      let all = feats ++ lib.optional profiling "PROFILING"; in
      lib.optionalString (all != []) ''export FEATURES="${lib.concatStringsSep " " all}"'';

    # Genuinely special environments, kept hand-written rather than in the matrix.
    # Each compute env is parameterized by profiling so it gets a -prof permutation;
    # emscripten (wasm) and codegen (build-only) have no meaningful profiling run.

    mkHeterogeneityCpu = { profiling ? false }: (pkgs.mkShell.override { stdenv = pkgs.stdenvNoCC; }) {
      name = "openlb-heterogeneity-cpu${lib.optionalString profiling "-prof"}";
      IN_OLB_NIX_SHELL = "1";
      buildInputs = baseTools ++ (with pkgs; [ gcc13 openmpi ]);
      shellHook = ''
        export ARCH="${marchExpr "tigerlake"}"
        export CXX=mpic++
        export CC=gcc

        export PLATFORMS="CPU_SISD CPU_SIMD"

        export CXXFLAGS="-O''${OLB_CXX_OPT_LEVEL:-3} -Wall -march=$ARCH -mtune=$ARCH -std=c++20"

        export PARALLEL_MODE=HYBRID
        export OMPFLAGS="-fopenmp"
        ${featuresLine [] profiling}
      '';
    };

    mkHeterogeneity = { profiling ? false }: (pkgs.mkShell.override { stdenv = pkgs.stdenvNoCC; }) {
      name = "openlb-heterogeneity${lib.optionalString profiling "-prof"}";
      IN_OLB_NIX_SHELL = "1";
      buildInputs = baseTools ++ cudaTools ++ (with pkgs; [ gcc13 openmpi-cuda ]);
      shellHook = ''
        export ARCH="${marchExpr "tigerlake"}"
        export CXX=mpic++
        export CC=gcc

        export PLATFORMS="CPU_SISD CPU_SIMD GPU_CUDA"

        export CXXFLAGS="-O3 -Wall -march=$ARCH -mtune=$ARCH -std=c++20"

        export PARALLEL_MODE=HYBRID
        export OMPFLAGS="-fopenmp"

        export CUDA_CXX=nvcc
        export CUDA_CXXFLAGS="-O3 -std=c++20 --forward-unknown-to-host-compiler"
        export CUDA_LDFLAGS="-L/run/opengl-driver/lib -L${cudaStubLib}"
        ${cudaArchLine}

        export FLOATING_POINT_TYPE=float
        ${featuresLine [] profiling}

        export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/run/opengl-driver/lib
      '';
    };

    specials = {
      heterogeneity-cpu = mkHeterogeneityCpu { };
      heterogeneity-cpu-prof = mkHeterogeneityCpu { profiling = true; };

      heterogeneity = mkHeterogeneity { };
      heterogeneity-prof = mkHeterogeneity { profiling = true; };

      emscripten = (pkgs.mkShell.override { stdenv = pkgs.stdenvNoCC; }) {
        name = "openlb-emscripten";
        IN_OLB_NIX_SHELL = "1";
        buildInputs = baseTools ++ (with pkgs; [ gcc13 emscripten ]);
        shellHook = ''
          export CXX=em++
          export CC=emcc

          export CXXFLAGS="-O3 -Wall -std=c++20 -Wno-missing-template-arg-list-after-template-kw"
          export LDFLAGS="-fPIC -g0 -s ALLOW_MEMORY_GROWTH -s MAXIMUM_MEMORY=4294967296 -s ENVIRONMENT=web -s EXPORTED_RUNTIME_METHODS=ccall,cwrap,UTF8ToString -fwasm-exceptions -lembind -s WASM=1"

          export PARALLEL_MODE=NONE
          export PLATFORMS="CPU_SISD"

          export FEATURES="EMSCRIPTEN"
        '';
      };

      # AMD GPU (HIP/ROCm) build environment. hipcc as sole driver compiling both
      # the CPU_SISD host code and the GPU_HIP device code (single-compiler mode,
      # like config/gpu_only_amd.mk). Targets an AMD gfx arch for compile coverage;
      # override the default via OLB_HIP_ARCH=<gfxXXX>. nixpkgs ROCm is AMD-only;
      # for HIP targeting NVIDIA use the hip-cuda shell below.
      hip = (pkgs.mkShell.override { stdenv = pkgs.stdenvNoCC; }) {
        name = "openlb-hip";
        IN_OLB_NIX_SHELL = "1";
        buildInputs = baseTools ++ (with pkgs; [
          binutils
          rocmPackages.clr rocmPackages.rocminfo
          rocmPackages.rocthrust rocmPackages.rocprim
        ]);
        shellHook = ''
          export CXX=hipcc
          export CC=hipcc

          export ROCM_PATH=${pkgs.rocmPackages.clr}
          export HIP_PLATFORM=amd
          export HIP_ARCH="''${OLB_HIP_ARCH:-gfx90a}"

          export PLATFORMS="CPU_SISD GPU_HIP"
          # --offload-arch and the host-only dependency-generation flag are added
          # by rules.mk based on HIP_PLATFORM/HIP_ARCH; only the rocThrust/rocPrim
          # include paths need to be injected here.
          export CXXFLAGS="-O3 -std=c++20 -I${pkgs.rocmPackages.rocthrust}/include -I${pkgs.rocmPackages.rocprim}/include"

          export PARALLEL_MODE=NONE
          export FLOATING_POINT_TYPE=float
        '';
      };

      # HIP targeting NVIDIA (HIP_PLATFORM=nvidia). Real ROCm hipcc forwards to
      # nvcc here, but nixpkgs' hipcc is hardwired to its AMD clang, so we drive
      # nvcc directly (which is what hipcc-nvidia does internally): the HIP API is
      # mapped onto CUDA via clr's base HIP headers plus the nvidia_detail headers
      # from ROCm/hipother. Everything is sourced from nixpkgs (clr headers,
      # cudaPackages for nvcc/runtime), so this is reproducible and runnable on an
      # NVIDIA GPU. The --generate-code / -x cu / CUDA-lib flags come from
      # rules.mk's HIP_PLATFORM=nvidia branch (all native nvcc options). Mirrors
      # config/gpu_only_hipcc_nvidia.mk.
      hip-cuda = pkgs.mkShell {
        name = "openlb-hip-cuda";
        IN_OLB_NIX_SHELL = "1";
        buildInputs = baseTools ++ cudaTools ++ (with pkgs; [
          binutils
          gcc13
          rocmPackages.clr
        ]);
        shellHook = ''
          export CXX=nvcc
          export CC=gcc

          export ROCM_PATH=${pkgs.rocmPackages.clr}
          export CUDA_PATH=${pkgs.cudaPackages.cudatoolkit}
          export HIP_PLATFORM=nvidia
          ${cudaArchLine}

          export PLATFORMS="CPU_SISD GPU_HIP"
          # Select the HIP NVIDIA platform and supply the CUDA runtime headers plus
          # the nvidia_detail HIP->CUDA mapping headers (absent from nixpkgs clr).
          # rules.mk adds the CUDA codegen/link flags for HIP_PLATFORM=nvidia.
          export CXXFLAGS="-O3 -std=c++20 -D__HIP_PLATFORM_NVIDIA__ -D__HIP_PLATFORM_NVCC__ -I${pkgs.cudaPackages.cudatoolkit}/include -I${hipother}/hipnv/include"
          export CUDA_LDFLAGS="-L/run/opengl-driver/lib -L${cudaStubLib}"

          export PARALLEL_MODE=NONE
          export FLOATING_POINT_TYPE=float

          # nixpkgs links libcuda against the CUDA driver *stub*; the real driver
          # lives in /run/opengl-driver/lib. Put it on the runtime search path so
          # HIP-on-CUDA binaries hit the real driver (else hipGetDeviceCount fails
          # with HIP error 999). Mirrors the native .#cuda shell — no per-job
          # LD_LIBRARY_PATH export needed.
          export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/run/opengl-driver/lib"
        '';
      };

      # Code generation / CSE optimization (needs mako + sympy).
      codegen = pkgs.mkShell {
        name = "openlb-codegen";
        IN_OLB_NIX_SHELL = "1";
        buildInputs = with pkgs; [
          gnumake
          gcc13
          gdb
          git
          (python3.withPackages (python-pkgs: with python-pkgs; [
            mako
            sympy
          ]))
        ];
        shellHook = ''
          export CXX=g++
          export CC=gcc

          export CXXFLAGS="-O0 -Wall -std=c++20"

          export PARALLEL_MODE=NONE

          export PLATFORMS="CPU_SISD"

          export FEATURES="INSPECT_DYNAMICS"
        '';
      };

      docs = pkgs.mkShell {
        name = "openlb-docs";
        buildInputs = with pkgs; let
          custom-texlive = pkgs.texlive.combine {
            inherit (pkgs.texlive)
              scheme-small xetex standalone pgfplots biblatex framed tcolorbox environ
              tikzfill pdfcol fontawesome5 algorithms algorithmicx placeins spath3
              siunitx algpseudocodex fifo-stack varwidth tabto-ltx totcount tikzmark
              tocloft
              # required by doc/include/latex_includes.tex
              scrhack setspaceenhanced lscapeenhanced floatbytocbasic;
          };

        in [
          quarto
          biber
          custom-texlive
          pdf2svg
          librsvg
          # the API reference is published next to the guide
          doxygen
          doxygen-awesome-css
          # doxygen_header.py, example_table.py
          python3
        ];

        shellHook = ''
          export DOXYGEN_AWESOME_CSS=${pkgs.doxygen-awesome-css}/share/doxygen-awesome-css
        '';
      };

      release-summary = pkgs.mkShell {
        name = "openlb-release-summary";
        buildInputs = with pkgs; [
          diff2html-cli
          curl
          git
          gnutar
        ];
      };
    };

  in {
    packages.${system} = {
      doxygen = pkgs.stdenvNoCC.mkDerivation rec {
        name = "openlb-doxygen";
        src = ./.;
        buildInputs = with pkgs; [
          gnumake
          doxygen
          doxygen-awesome-css
        ];
        buildPhase = ''
          make doxygen DOXYGEN_AWESOME_CSS=${pkgs.doxygen-awesome-css}/share/doxygen-awesome-css
        '';
        installPhase = ''
          mkdir -p $out
          cp -r build/doxygen/html $out/
        '';
      };
    };

    devShells.${system} = generated // specials // {
      default = generated.gcc;
      # the docs environment was called quarto before it also built the API reference
      quarto = specials.docs;
    };
  };
}

---
name: mpi-domain-decomposition
description: Implement or review MPI domain decomposition in this repository, especially slab decomposition, halo exchange, reductions, and gathered diagnostics in `3DRBOpenaccMpi.F90` or derived MPI solvers. Use when adding ranks, `MPI_Sendrecv`, halo layers, `MPI_Allreduce`, `MPI_Gather`, global or local indexing, or rank-specific output.
---

# MPI Domain Decomposition

Use this skill when a change affects how the solver is split across ranks or how rank-local results are assembled back into global results.

## Repository baseline

Use `3DRBOpenaccMpi.F90` as the source of truth.

The current MPI design is:

- 1D `z`-slab decomposition
- owned cells on `1:nz`
- halo planes on `0` and `nz+1`
- global range tracked by `zStartGlobal` and `zEndGlobal`
- neighbor ranks tracked by `mpiLeft` and `mpiRight`

Do not silently change this decomposition model. Call out any move away from 1D `z`-slabs as an architectural departure.

## Default workflow

### 1. Define the local subdomain correctly

Update the split logic together:

- `nzGlobal`, `nz`, `zStartGlobal`, `zEndGlobal`
- first-rank and last-rank flags
- periodic-neighbor handling
- root-rank assumptions for reporting

Keep local sizes strictly positive and guard impossible rank counts early.

### 2. Keep halo ownership unambiguous

Treat:

- `1:nz` as owned interior planes
- `0` as the lower halo
- `nz+1` as the upper halo

When porting a new field into MPI, update allocation, initialization, reload, halo exchange, and post-processing together.

### 3. Keep communication pairs symmetric

Match send and receive logic in one place:

- same buffer length on both sides
- distinct tags for distinct data families
- `MPI_PROC_NULL` branches for non-periodic outer ranks
- identical pack and unpack traversal order

For the current repository, mirror the style of `exchange_f_post_halo_mpi` and `exchange_g_post_halo_mpi`.

### 4. Keep global reductions correct

Use:

- `MPI_Allreduce` for global sums, residuals, `Nu`, and `Re`
- `MPI_Gather` or `MPI_Gatherv` when root must choose a maximum or reconstruct ordered data

Do not trust a global diagnostic until the reduction path is verified.

### 5. Keep local and global coordinates distinct

When computing extrema or mid-plane quantities:

- use local coordinates only for rank-local search
- convert or compare against global coordinates before reporting
- handle the case where the target global plane is not owned by a given rank

Mirror the style of `calc_umid_max_common3d`, `calc_vmid_max_common3d`, and `calc_wmid_max_common3d`.

## MPI plus OpenACC rule

If arrays live on the device:

- update the boundary planes back to host before MPI communication
- perform MPI on host buffers unless the code explicitly supports device-aware MPI
- update the modified halo data back to the device before the next GPU kernel

Use `openacc-porting` together with this skill when both concerns are active.

## Common failure modes

- Mix local `nz` with global `nzGlobal`.
- Forget to update both flow and thermal halo exchanges.
- Reuse a tag and accidentally cross-wire two message families.
- Write root-only files from all ranks.
- Validate local `Nu` or residuals but forget the final `MPI_Allreduce`.

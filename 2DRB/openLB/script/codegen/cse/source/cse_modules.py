# This file is part of the OpenLB library
#
# Copyright (C) 2024 Shota Ito, Adrian Kummerlaender
# E-mail contact: info@openlb.net
# The most recent release of OpenLB can be downloaded at
# <http://www.openlb.net/>
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public
# License along with this program; if not, write to the Free
# Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
# Boston, MA  02110-1301, USA.

from mako.template import Template
from pathlib import Path
from metadata import getHardwareMetadata, MakefileEnvironment, getCommitHash
import subprocess, os
import json

def compileDynamics(dynamics_type, template, dynamics_hash, work_dir, precision="float"):
    work_path = Path(work_dir)
    template_extraction = Template(filename=template)
    
    executable = dynamics_hash
    output_cpp = work_path / f"{dynamics_hash}.cpp"
    output_template = work_path / f"{dynamics_hash}.out"

    base_dir = Path.cwd() 
    tmp_dir = base_dir / "tmp_results"
    json_path = tmp_dir / f"result_{dynamics_hash}.json"
    
    compileFlags= [
        "CXXFLAGS", 
        "CXX", 
        "CC", 
        "PARALLEL_MODE", 
        "PLATFORMS", 
        "OMPFLAGS", 
        "FEATURES"
    ]
    
    m = MakefileEnvironment(work_dir)
    compile_flags_json = m.get(compileFlags)
     
    hardware_metadata = getHardwareMetadata(compile_flags_json["PLATFORMS"])
    commit_hash = getCommitHash()

    # fill mako template
    with open(output_cpp, 'w') as f:
        f.write(
                template_extraction.render(
                    dynamics=dynamics_type, 
                    output=output_template, 
                    env=json.dumps(compile_flags_json, indent=8),
                    system=json.dumps(hardware_metadata, indent=8),
                    commitHash=commit_hash,
                    precision=precision,
                    jsonPath=str(json_path.absolute())
                )
        )
    
    # compile and execute
    subprocess.call("make EXAMPLE="+executable+" > /dev/null", shell=True, cwd=work_dir)
  
def runProfileDynamics(dynamics_hash, work_dir, size=None, runs=None):
    
    executable = dynamics_hash;

    executable_args = ""
    executable_args += f"--nCells {size}" if size != None else ""
    executable_args += f" --runs {runs}" if runs   != None else ""
    
    print("Start Profiling")
    subprocess.call(f"rocprof -i ../../counters.txt --timestamp on -o roofline_data.csv ./{executable} {executable_args}", shell=True, cwd=work_dir)
    

def runDynamics(dynamics_hash, work_dir, size=None, runs=None):

    executable = f"./{dynamics_hash}"
    args = []
    if size: args.extend(["--nCells", str(size)])
    if runs: args.extend(["--runs", str(runs)])
    
    subprocess.run([executable] + args, cwd=work_dir, check=True)

def benchmarkDynamics(dynamics_type, template, dynamics_hash, work_dir, size=None, runs=None, precision="float"):

    template_extraction = Template(filename=template)
    base_dir = Path.cwd() 
    tmp_dir = base_dir / "tmp_results"
    json_path = tmp_dir / f"result_{dynamics_hash}.json"
    
    # create subdirectories to prevent confusion by linker
    executable = dynamics_hash;

    # specify output for the cpp mako template
    output_cpp = work_dir + dynamics_hash + ".cpp"
    output_template = work_dir + dynamics_hash + ".out"

    compileFlags= [
        "CXXFLAGS", 
        "CXX", 
        "CC", 
        "PARALLEL_MODE", 
        "PLATFORMS", 
        "OMPFLAGS", 
        "FEATURES"
    ]
    
    m = MakefileEnvironment(work_dir)
    compile_flags_json = m.get(compileFlags)
     
    hardware_metadata = getHardwareMetadata(compile_flags_json["PLATFORMS"])
    commit_hash = getCommitHash()

    # fill mako template
    with open(output_cpp, 'w') as f:
        f.write(
                template_extraction.render(
                    dynamics=dynamics_type, 
                    output=output_template, 
                    env=json.dumps(compile_flags_json, indent=8),
                    system=json.dumps(hardware_metadata, indent=8),
                    commitHash=commit_hash,
                    precision=precision,
                    jsonPath=str(json_path.absolute())
                )
        )
    
    # compile and execute
    subprocess.call("make EXAMPLE="+executable+" > /dev/null", shell=True, cwd=work_dir)
    
    executable_args = ""
    executable_args += f"--nCells {size}" if size != None else ""
    executable_args += f" --runs {runs}" if runs   != None else ""
 
    subprocess.call(f"./{executable} {executable_args}", shell=True, cwd=work_dir)

def extractDynamicsExpressions(dynamics_type, template, input_f, output_p, output_o, disable_cse=False, precision="float", target="Generic"):
    # Fill Mako template of C++ program to retrieve expression tree
    template_extraction = Template(filename=template)
    base_dir = Path.cwd()

    # create subdirectories to prevent confusion by linker
    executable = Path(input_f).stem;

    # specify output for the cpp mako template
    output_cpp = output_p + Path(input_f).stem + ".cpp"
    output_template = output_o + Path(input_f).stem + ".out"

    # fill mako template
    with open(output_cpp, 'w') as f:
        f.write(
                template_extraction.render(
                    dynamics=dynamics_type,
                    output=output_template,
                    disable_cse=disable_cse,
                    precision=precision,
                    target=target
                )
        )
    
    # compile and execute
    subprocess.call("make EXAMPLE="+executable+" > /dev/null", shell=True, cwd=output_p)
    
    executable_args = ""
    subprocess.call(f"./{executable} {executable_args}", shell=True, cwd=output_p)

def extractOperatorExpressions(operator_with_descriptor_type, template, input_f, output_p, output_o): # Fill Mako template of C++ program to retrieve expression tree
    template_extraction = Template(filename=template)
   
    # create subdirectories to prevent confusion by linker
    executable = Path(input_f).stem;

    # specify output for the cpp mako template
    output_cpp = output_p + Path(input_f).stem + ".cpp"
    output_template = output_o + Path(input_f).stem + ".out"

    # split operator and descriptor
    operator_type, descriptor_type = operator_with_descriptor_type.split(';',1)

    # fill mako template
    with open(output_cpp, 'w') as f:
        f.write(template_extraction.render(operator=operator_type, descriptor=descriptor_type, output=output_template))
    
    # compile and execute
    subprocess.call("make EXAMPLE="+executable+";", shell=True, cwd=output_p)
    subprocess.call("./"+executable, shell=True, cwd=output_p)

def executeCSE(expressions, template, result, deterministic=False, target="Generic"):
    # Generate cse-optimized C++ code from tree and fill template
    template_cse = Template(filename=template)
    with open(result, 'w') as f:
        f.write(template_cse.render(filename=expressions, deterministic=deterministic, target=target))

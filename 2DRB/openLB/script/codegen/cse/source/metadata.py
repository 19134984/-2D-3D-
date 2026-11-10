import subprocess, os
import re

class MakefileEnvironment():
    def __init__(self, target_dir):
        self.dryRun = subprocess.run(
            ['make', '-p', '--dry-run'],
            capture_output=True,
            text=True,
            cwd=target_dir
        ).stdout

    def get(self, variables):
        flags={}
        for line in self.dryRun.splitlines():
            if not line.startswith('\t') and '=' in line:
                left, right = line.split('=', 1)
                
                var_name = left.strip(' \t:+?!')
                if var_name in variables:
                    flags[var_name] = right.split()
        return flags

def getHardwareMetadata(platforms):
    def getCPUModel():
        try:
        # Reads the model name directly from the system's CPU info
            with open("/proc/cpuinfo", "r") as f:
                for line in f:
                    if "model name" in line:
                        return re.sub( ".*model name.*:", "", line, 1).strip()
        except Exception:
            return "Unknown CPU"
    def getNvidiaModel():
        try:
            # --query-gpu=name returns just the model (e.g., "NVIDIA A100-SXM4-80GB")
            result = subprocess.run(
                ["nvidia-smi", "--query-gpu=name", "--format=csv,noheader"],
                capture_output=True, text=True, check=True
            )
            return result.stdout.strip()
        except:
            return None
    def getAMDModel():
        try:
            result = subprocess.run(
                ["rocm-smi", "--showproductname"],
                capture_output=True, text=True, check=True
            )
            for line in result.stdout.splitlines():
                if "Card Series" in line or "Product Name" in line:
                    return line.split(":")[-1].strip()
        except:
            return None
    def getIntelModel():
        try:
            result = subprocess.run(
                ["sycl-ls"],
                capture_output=True, text=True, check=True
            )
            for line in result.stdout.splitlines():
                if "gpu" in line.lower():
                    return line.split("]", 1)[-1].strip()
        except:
            try:
                result = subprocess.run(
                    ["xpu-smi", "discovery", "--list"],
                    capture_output=True, text=True, check=True
                )
                for line in result.stdout.splitlines():
                    if "Model Name" in line:
                        return line.split(":")[-1].strip()
            except:
                return None
        return None 
    gpu = "None"
    
    if "GPU_CUDA" in platforms:
        gpu = getNvidiaModel()
    if "GPU_HIP" in platforms:
        gpu = getAMDModel()
    if "GPU_SYCL" in platforms:
        gpu = getIntelModel() or gpu
    cpu = getCPUModel()
    
    return {
        "cpu": cpu,
        "gpu": gpu,
        "os": subprocess.run(["uname", "-sr"], capture_output=True, text=True).stdout.strip()
    }



def getCommitHash(short=True):
    cmd = ["git", "rev-parse", "HEAD"]
    if short:
        cmd.insert(2, "--short")

    return subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            ).stdout.strip()


# pySCA UX Implementation Guide

This guide provides detailed implementation specifications for the UX improvements outlined in UX_DESIGN.md.

---

## 1. Progress Indicators

### 1.1 Requirements
- Use `tqdm` library for progress bars
- Show progress for operations > 5 seconds
- Display nested progress for multi-step operations
- Provide time estimates for long operations

### 1.2 Implementation Pattern

```python
# pysca/progress.py
from tqdm import tqdm
import time
from typing import Optional, Callable

class ProgressTracker:
    """Context manager for tracking progress of operations."""
    
    def __init__(self, description: str, total: Optional[int] = None, 
                 unit: str = "items", disable: bool = False):
        self.description = description
        self.total = total
        self.unit = unit
        self.disable = disable
        self.pbar = None
    
    def __enter__(self):
        self.pbar = tqdm(
            total=self.total,
            desc=self.description,
            unit=self.unit,
            disable=self.disable,
            bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]'
        )
        return self
    
    def __exit__(self, *args):
        if self.pbar:
            self.pbar.close()
    
    def update(self, n: int = 1):
        if self.pbar:
            self.pbar.update(n)
    
    def set_description(self, desc: str):
        if self.pbar:
            self.pbar.set_description(desc)

# Usage example in scaProcessMSA_py3_big.py
from pysca.progress import ProgressTracker

def process_alignment(sequences, logger, quiet=False):
    with ProgressTracker("Processing sequences", total=len(sequences), disable=quiet) as pbar:
        for i, seq in enumerate(sequences):
            # ... process sequence ...
            pbar.update(1)
            if i % 1000 == 0:
                pbar.set_description(f"Processing sequences ({i}/{len(sequences)})")
```

### 1.3 Integration Points

**scaProcessMSA_py3_big.py:**
- MMseqs2 preclustering progress (use subprocess output parsing)
- Sequence filtering progress
- Weight computation progress

**scaCore_py3.py:**
- Matrix computation progress (use chunks)
- Randomization trials progress
- Sequence correlation computation progress

---

## 2. Preset Configurations

### 2.1 Preset Definitions

```python
# pysca/presets.py
from dataclasses import dataclass
from typing import Optional

@dataclass
class PreprocessingPreset:
    """Preset configuration for preprocessing step."""
    max_gap_pos: float
    max_gap_seq: float
    min_SID: float
    max_SID: float
    initial_trim_gap: float
    precluster: Optional[bool]  # None = auto
    cluster_id: float
    cluster_coverage: float

@dataclass
class CorePreset:
    """Preset configuration for core SCA step."""
    lbda: float
    norm: str
    Ntrials: int
    do_seqcorr: bool
    do_sector_id: bool
    sector_cutoff: float
    float32: bool

@dataclass
class WorkflowPreset:
    """Complete workflow preset."""
    name: str
    description: str
    preprocessing: PreprocessingPreset
    core: CorePreset

# Preset definitions
PRESETS = {
    "quick": WorkflowPreset(
        name="quick",
        description="Fast analysis with minimal parameters",
        preprocessing=PreprocessingPreset(
            max_gap_pos=0.3,  # More lenient
            max_gap_seq=0.3,
            min_SID=0.15,
            max_SID=0.9,
            initial_trim_gap=0.9,
            precluster=True,  # Always precluster for speed
            cluster_id=0.8,
            cluster_coverage=0.75
        ),
        core=CorePreset(
            lbda=0.03,
            norm="frob",
            Ntrials=5,  # Fewer trials
            do_seqcorr=False,  # Skip for speed
            do_sector_id=False,  # Skip for speed
            sector_cutoff=0.95,
            float32=True  # Use float32 for speed
        )
    ),
    
    "standard": WorkflowPreset(
        name="standard",
        description="Balanced analysis with recommended settings",
        preprocessing=PreprocessingPreset(
            max_gap_pos=0.2,
            max_gap_seq=0.2,
            min_SID=0.2,
            max_SID=0.8,
            initial_trim_gap=0.8,
            precluster=None,  # Auto-detect
            cluster_id=0.85,
            cluster_coverage=0.8
        ),
        core=CorePreset(
            lbda=0.03,
            norm="frob",
            Ntrials=10,
            do_seqcorr=True,
            do_sector_id=True,
            sector_cutoff=0.95,
            float32=False
        )
    ),
    
    "detailed": WorkflowPreset(
        name="detailed",
        description="Comprehensive analysis with strict filtering",
        preprocessing=PreprocessingPreset(
            max_gap_pos=0.15,  # Stricter
            max_gap_seq=0.15,
            min_SID=0.25,
            max_SID=0.75,
            initial_trim_gap=0.7,
            precluster=False,  # No preclustering for accuracy
            cluster_id=0.9,
            cluster_coverage=0.85
        ),
        core=CorePreset(
            lbda=0.03,
            norm="frob",
            Ntrials=20,  # More trials
            do_seqcorr=True,
            do_sector_id=True,
            sector_cutoff=0.95,
            float32=False
        )
    )
}

def get_preset(name: str) -> WorkflowPreset:
    """Get preset by name."""
    if name not in PRESETS:
        raise ValueError(f"Unknown preset: {name}. Available: {list(PRESETS.keys())}")
    return PRESETS[name]

def list_presets() -> list[str]:
    """List available presets."""
    return list(PRESETS.keys())
```

### 2.2 Integration with Argument Parsers

```python
# In scaProcessMSA_py3_big.py main()
def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser()
    # ... existing arguments ...
    
    # Add preset argument
    p.add_argument("--preset", choices=["quick", "standard", "detailed"],
                   help="Use preset configuration (overrides other parameters)")
    
    args = p.parse_args(argv)
    
    # Apply preset if specified
    if args.preset:
        from pysca.presets import get_preset
        preset = get_preset(args.preset)
        preproc = preset.preprocessing
        
        # Override args with preset values (unless explicitly set)
        if not hasattr(args, '_explicit_parameters'):
            args.parameters = [
                preproc.max_gap_pos, preproc.max_gap_seq,
                preproc.min_SID, preproc.max_SID
            ]
        if not hasattr(args, '_explicit_precluster'):
            args.precluster = preproc.precluster
        if not hasattr(args, '_explicit_cluster_id'):
            args.cluster_id = preproc.cluster_id
        # ... etc for other parameters ...
    
    # Continue with normal processing...
```

---

## 3. Configuration Files

### 3.1 Configuration File Format (YAML)

```yaml
# Example: sca_standard.yaml
name: "standard_analysis"
description: "Standard SCA analysis with balanced parameters"

preprocessing:
  reference:
    pdb: "1XYZ"
    chain: "A"
  
  filtering:
    max_gap_pos: 0.2
    max_gap_seq: 0.2
    min_SID: 0.2
    max_SID: 0.8
    initial_trim_gap: 0.8
  
  preclustering:
    enabled: null  # null = auto-detect
    cluster_id: 0.85
    cluster_coverage: 0.8

core:
  regularization: 0.03
  norm: "frob"
  trials: 10
  sequence_correlations: true
  sector_identification: true
  sector_cutoff: 0.95
  kpos: 0  # auto
  float32: false

output:
  matlab: false
  save_msa_numeric: false
```

### 3.2 Configuration File Parser

```python
# pysca/config.py
import yaml
from pathlib import Path
from typing import Dict, Any, Optional
from dataclasses import dataclass, asdict

@dataclass
class Config:
    """Configuration container."""
    name: str
    description: str
    preprocessing: Dict[str, Any]
    core: Dict[str, Any]
    output: Dict[str, Any]
    
    @classmethod
    def from_yaml(cls, path: Path) -> 'Config':
        """Load configuration from YAML file."""
        with open(path, 'r') as f:
            data = yaml.safe_load(f)
        return cls(**data)
    
    def to_yaml(self, path: Path) -> None:
        """Save configuration to YAML file."""
        with open(path, 'w') as f:
            yaml.dump(asdict(self), f, default_flow_style=False, sort_keys=False)
    
    @classmethod
    def from_args(cls, args) -> 'Config':
        """Create config from command-line arguments."""
        # Convert argparse args to config structure
        # This allows saving current run as config
        pass

def load_config(path: Path) -> Config:
    """Load configuration file."""
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")
    return Config.from_yaml(path)

def save_config(config: Config, path: Path) -> None:
    """Save configuration file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    config.to_yaml(path)

def validate_config(config: Config) -> list[str]:
    """Validate configuration and return list of warnings/errors."""
    errors = []
    
    # Validate preprocessing
    if 'preprocessing' in config.__dict__:
        preproc = config.preprocessing
        if 'filtering' in preproc:
            filt = preproc['filtering']
            if filt.get('max_gap_pos', 0) > 1.0 or filt.get('max_gap_pos', 0) < 0:
                errors.append("max_gap_pos must be between 0 and 1")
            # ... more validation ...
    
    return errors
```

### 3.3 Integration

```python
# In main() functions
p.add_argument("--config", type=Path, help="Load configuration from YAML file")
p.add_argument("--save-config", type=Path, help="Save current parameters as config file")

args = p.parse_args(argv)

# Load config if provided
if args.config:
    from pysca.config import load_config, validate_config
    config = load_config(args.config)
    errors = validate_config(config)
    if errors:
        logger.error("Config validation errors:")
        for err in errors:
            logger.error(f"  - {err}")
        return 1
    
    # Apply config to args
    # ... merge config values into args ...
```

---

## 4. Enhanced Error Messages

### 4.1 Error Message Framework

```python
# pysca/errors.py
from enum import Enum
from typing import Optional, List

class ErrorType(Enum):
    FILE_NOT_FOUND = "file_not_found"
    INVALID_FORMAT = "invalid_format"
    MISSING_DEPENDENCY = "missing_dependency"
    INVALID_PARAMETER = "invalid_parameter"
    MEMORY_ERROR = "memory_error"
    WORKFLOW_ERROR = "workflow_error"

class SCAAwareError(Exception):
    """Base exception with user-friendly error messages."""
    
    def __init__(self, message: str, error_type: ErrorType,
                 suggestions: Optional[List[str]] = None,
                 documentation_link: Optional[str] = None):
        self.message = message
        self.error_type = error_type
        self.suggestions = suggestions or []
        self.documentation_link = documentation_link
        super().__init__(self.format_message())
    
    def format_message(self) -> str:
        """Format error message with suggestions."""
        lines = [f"❌ Error: {self.message}"]
        
        if self.suggestions:
            lines.append("\nPossible solutions:")
            for i, suggestion in enumerate(self.suggestions, 1):
                lines.append(f"  {i}. {suggestion}")
        
        if self.documentation_link:
            lines.append(f"\nSee: {self.documentation_link}")
        
        return "\n".join(lines)

# Specific error classes
class PDBNotFoundError(SCAAwareError):
    def __init__(self, pdb_id: str):
        super().__init__(
            f"PDB file not found: {pdb_id}",
            ErrorType.FILE_NOT_FOUND,
            suggestions=[
                f"Verify PDB ID at https://www.rcsb.org/structure/{pdb_id}",
                "Check internet connection (PDB files are downloaded automatically)",
                f"Provide local path: --pdb /path/to/{pdb_id}.pdb",
                "Use --refseq to provide a reference sequence file instead"
            ],
            documentation_link="https://ranganathanlab.gitlab.io/pySCA/usage.html#reference-sequence"
        )

class LargeMSAWarning(SCAAwareError):
    def __init__(self, n_sequences: int):
        super().__init__(
            f"Large MSA detected ({n_sequences:,} sequences)",
            ErrorType.WORKFLOW_ERROR,
            suggestions=[
                "Consider using --precluster for faster processing",
                f"Expected processing time without preclustering: {n_sequences // 1000} hours",
                "Use --preset standard for automatic optimization",
                "See LARGE_MSA_OPTIMIZATIONS.md for details"
            ]
        )
```

### 4.2 Usage in Code

```python
# In scaProcessMSA_py3_big.py
from pysca.errors import PDBNotFoundError, LargeMSAWarning

def load_pdb(pdb_id: str, chain: str) -> Structure:
    try:
        # ... attempt to load PDB ...
    except FileNotFoundError:
        raise PDBNotFoundError(pdb_id) from None

def analyze_msa_size(sequences: List[str], logger) -> None:
    n_seq = len(sequences)
    if n_seq > 100000:
        warning = LargeMSAWarning(n_seq)
        logger.warning(str(warning))
```

---

## 5. Input Validation

### 5.1 Pre-flight Validation

```python
# pysca/validation.py
from pathlib import Path
from typing import List, Tuple, Optional
import sys

class ValidationResult:
    """Result of validation check."""
    def __init__(self, valid: bool, errors: List[str], warnings: List[str]):
        self.valid = valid
        self.errors = errors
        self.warnings = warnings
    
    def print(self) -> None:
        """Print validation results."""
        if self.errors:
            print("❌ Validation Errors:", file=sys.stderr)
            for err in self.errors:
                print(f"  - {err}", file=sys.stderr)
        
        if self.warnings:
            print("⚠️  Validation Warnings:", file=sys.stderr)
            for warn in self.warnings:
                print(f"  - {warn}", file=sys.stderr)
        
        if self.valid and not self.warnings:
            print("✓ Validation passed")

def validate_alignment_file(path: Path) -> ValidationResult:
    """Validate alignment file before processing."""
    errors = []
    warnings = []
    
    # Check file exists
    if not path.exists():
        errors.append(f"Alignment file not found: {path}")
        return ValidationResult(False, errors, warnings)
    
    # Check file is readable
    if not os.access(path, os.R_OK):
        errors.append(f"Cannot read alignment file: {path}")
        return ValidationResult(False, errors, warnings)
    
    # Check file size
    size_mb = path.stat().st_size / (1024 * 1024)
    if size_mb > 1000:
        warnings.append(f"Large alignment file ({size_mb:.1f} MB) - processing may be slow")
    
    # Try to parse file
    try:
        from Bio import AlignIO
        alignment = AlignIO.read(str(path), "fasta")
        n_seq = len(alignment)
        n_pos = alignment.get_alignment_length()
        
        if n_seq < 10:
            errors.append(f"Too few sequences ({n_seq}) - need at least 10")
        
        if n_pos < 20:
            errors.append(f"Alignment too short ({n_pos} positions) - need at least 20")
        
        if n_seq > 500000:
            warnings.append(f"Very large MSA ({n_seq:,} sequences) - consider preclustering")
    
    except Exception as e:
        errors.append(f"Cannot parse alignment file: {e}")
    
    return ValidationResult(len(errors) == 0, errors, warnings)

def validate_pdb(pdb_id: str, chain: str) -> ValidationResult:
    """Validate PDB structure."""
    errors = []
    warnings = []
    
    # Check if PDB ID format is valid
    if len(pdb_id) != 4:
        errors.append(f"Invalid PDB ID format: {pdb_id} (should be 4 characters)")
    
    # Try to fetch/access PDB
    # ... implementation ...
    
    return ValidationResult(len(errors) == 0, errors, warnings)

def validate_workflow_preprocessing(alignment_path: Path, pdb_id: Optional[str] = None) -> ValidationResult:
    """Comprehensive validation for preprocessing step."""
    errors = []
    warnings = []
    
    # Validate alignment
    align_result = validate_alignment_file(alignment_path)
    errors.extend(align_result.errors)
    warnings.extend(align_result.warnings)
    
    # Validate PDB if provided
    if pdb_id:
        pdb_result = validate_pdb(pdb_id, "A")
        errors.extend(pdb_result.errors)
        warnings.extend(pdb_result.warnings)
    
    # Check disk space
    # ... implementation ...
    
    # Check memory (rough estimate)
    # ... implementation ...
    
    return ValidationResult(len(errors) == 0, errors, warnings)
```

### 5.2 Integration

```python
# In main() functions
p.add_argument("--validate", action="store_true",
               help="Validate inputs and exit without processing")

args = p.parse_args(argv)

if args.validate:
    from pysca.validation import validate_workflow_preprocessing
    result = validate_workflow_preprocessing(
        Path(args.alignment),
        args.pdb if hasattr(args, 'pdb') else None
    )
    result.print()
    return 0 if result.valid else 1
```

---

## 6. Summary Output

### 6.1 Summary Generator

```python
# pysca/summary.py
from pathlib import Path
from datetime import datetime
from typing import Dict, Any, Optional
import json

def generate_summary(data: Dict[str, Any], output_dir: Path) -> Path:
    """Generate human-readable summary file."""
    summary_path = output_dir / "summary.txt"
    
    lines = [
        "=" * 60,
        "SCA Analysis Summary",
        "=" * 60,
        f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "",
    ]
    
    if 'input' in data:
        lines.extend([
            "## Input",
            f"  Alignment file: {data['input'].get('alignment', 'N/A')}",
            f"  Reference: {data['input'].get('reference', 'N/A')}",
            "",
        ])
    
    if 'preprocessing' in data:
        preproc = data['preprocessing']
        lines.extend([
            "## Preprocessing Results",
            f"  Input sequences: {preproc.get('input_sequences', 'N/A'):,}",
            f"  Final sequences: {preproc.get('final_sequences', 'N/A'):,}",
            f"  Positions: {preproc.get('positions', 'N/A')}",
            f"  Effective sequences: {preproc.get('effective_sequences', 'N/A'):.0f}",
            "",
        ])
    
    if 'core' in data:
        core = data['core']
        lines.extend([
            "## SCA Core Results",
            f"  Regularization (λ): {core.get('lbda', 'N/A')}",
            f"  Randomization trials: {core.get('trials', 'N/A')}",
            "",
        ])
    
    if 'sector' in data:
        sector = data['sector']
        lines.extend([
            "## Sector Identification",
            f"  Sectors identified: {sector.get('n_sectors', 'N/A')}",
            f"  Eigenmodes used: {sector.get('kpos', 'N/A')}",
            "",
        ])
    
    lines.extend([
        "## Output Files",
    ])
    for file_type, path in data.get('output_files', {}).items():
        lines.append(f"  {file_type}: {path}")
    
    lines.extend([
        "",
        "=" * 60,
    ])
    
    with open(summary_path, 'w') as f:
        f.write('\n'.join(lines))
    
    return summary_path
```

---

## 7. Unified Workflow Command

### 7.1 Workflow Runner

```python
# pysca/workflow.py
from pathlib import Path
from typing import Optional
import subprocess
import sys

def run_complete_workflow(
    alignment: Path,
    pdb_id: Optional[str] = None,
    chain: str = "A",
    preset: str = "standard",
    output_name: Optional[str] = None,
    **kwargs
) -> int:
    """Run complete SCA workflow in one command."""
    
    # Step 1: Preprocessing
    cmd1 = [
        "python", "-m", "pysca.scaProcessMSA_py3_big",
        str(alignment)
    ]
    
    if pdb_id:
        cmd1.extend(["-s", pdb_id, "--chainID", chain])
    
    if preset:
        cmd1.extend(["--preset", preset])
    
    if output_name:
        cmd1.extend(["--output", output_name])
    
    print("Step 1/2: Preprocessing alignment...")
    result1 = subprocess.run(cmd1)
    if result1.returncode != 0:
        return result1.returncode
    
    # Determine output database path
    if output_name:
        db_path = Path("Outputs") / f"{output_name}.db.gz"
    else:
        db_path = Path("Outputs") / f"{alignment.stem}.db.gz"
    
    # Step 2: Core analysis
    cmd2 = [
        "python", "-m", "pysca.scaCore_py3",
        str(db_path),
        "--preset", preset
    ]
    
    print("Step 2/2: Running SCA core calculations...")
    result2 = subprocess.run(cmd2)
    
    return result2.returncode

# CLI entry point
# scripts/sca-run
#!/usr/bin/env python3
"""Unified workflow command."""
import sys
from pathlib import Path
from pysca.workflow import run_complete_workflow

def main():
    # Parse arguments...
    result = run_complete_workflow(...)
    sys.exit(result)

if __name__ == "__main__":
    main()
```

---

## Implementation Checklist

- [ ] Add `tqdm` dependency to setup.py
- [ ] Create `pysca/progress.py` with ProgressTracker
- [ ] Create `pysca/presets.py` with preset definitions
- [ ] Create `pysca/config.py` for configuration file handling
- [ ] Create `pysca/errors.py` with enhanced error classes
- [ ] Create `pysca/validation.py` for input validation
- [ ] Create `pysca/summary.py` for summary generation
- [ ] Integrate progress tracking in scaProcessMSA_py3_big.py
- [ ] Integrate progress tracking in scaCore_py3.py
- [ ] Add preset support to argument parsers
- [ ] Add configuration file support
- [ ] Add validation mode (`--validate`)
- [ ] Create unified workflow command (`sca-run`)
- [ ] Add enhanced error messages throughout
- [ ] Update documentation with new features



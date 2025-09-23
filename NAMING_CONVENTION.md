## Result File Naming Convention

After running jobs, results will be saved with **distinct method prefixes** for easy identification:

### File Name Pattern:
- **Standard**: `standard_L8_lam0.1_s1.json`
- **Forced**: `forced_L8_lam0.1_s1.json` 
- **Dummy**: `dummy_L8_lam0.1_s1.json`

### Easy Result Collection:
```bash
# Count completed jobs by method
ls output/standard_*.json | wc -l    # Standard BinderSim results
ls output/forced_*.json | wc -l      # Forced BinderSim results  
ls output/dummy_*.json | wc -l       # Dummy site results

# Analyze specific method
python analyze.py output/standard_*.json    # Analyze all standard results
python analyze.py output/forced_*.json      # Analyze all forced results
python analyze.py output/dummy_*.json       # Analyze all dummy results
```

### Submission Commands:
- `./submit_standard.sh` → Creates `output/standard_*.json`
- `./submit_forced.sh` → Creates `output/forced_*.json`  
- `./submit_dummy.sh` → Creates `output/dummy_*.json`

This makes it super easy to find and analyze results from each method separately!
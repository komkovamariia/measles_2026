# Repository validation

Run

```bash
python scripts/validate_repository.py
```

The validator checks

- 112 sequencing libraries are represented
- 96 new libraries and 16 legacy libraries are present
- T0 is absent
- controls use CONTROL only
- child means age below 18 years
- adult means age 18 years and older
- final library QC threshold is 4400 productive TRB clonotypes
- final cohort contains 78 unique patients, 33 children and 45 adults
- cross sectional counts are 4 controls, 27 T1 and 8 T2 for children
- cross sectional counts are 13 controls, 29 T1 and 19 T2 for adults
- paired T1 and T2 counts are 6 children and 16 adults
- `new_111` is included as a T2 biological repeat of `new_106`
- the old age 6 reporting script is absent
- expanded clone analysis uses the 0.10 percent threshold
- the corrected legacy patient identities are preserved

The validator checks repository configuration and analysis contracts. It does not replace a successful analysis run on the cluster

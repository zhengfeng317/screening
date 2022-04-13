# screening

1. go to POT/; change VASP pp path in get-pot.sh then execute 
2. go to str/; change RE elements and natom
3. go to NM; change init-rlx/job.pbs ; change 'POT' path in the  0-sub-rlx.sh
- submit NM
- check via check-NM.py

4. go to FM; change "MAGMOM" in INCAR.* files (set metal moment=1; others =0)

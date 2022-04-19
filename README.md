# screening

1. go to POT/; change VASP pp path in get-pot.sh then execute 
2. go to str/; change RE elements and natom
3. go to NM; change init-rlx/job.pbs ; change 'POT' path in the  0-sub-rlx.sh
- submit NM
- check via check-NM.py

4. go to FM; change "MAGMOM" in INCAR.* files (set metal moment=1; others =0)
- submit FM via 1-sub-fm.sh 
- check via check-NM.py

5. Once finished, get results by results-FM-NM.py, then plot by plot-FM-NM.py.
6. Analyse FM-NM.png for either EP or SF calculations.
7. go to AFM; provide system for '1-sub-afm.sh'; change "MAGMOM" in INCAR.* files with AFM configurations
- submit AFM via 1-sub-afm.sh 
- analyze results via results-SF.py

8. go to EP; provide system for '2-sub-ep.sh'; 
- submit EP via 2-sub-ep.sh 
- check and resubmit unfinished jobs with check.sh/resub.sh (may need to modify)
- use process-phonon.sh to process the phonon results; process-omega.sh to compute omega and plot

9. go to CVX to compute convex hull for structures of interest; provide systems for "0-get-str.sh"
- download MP referece list via 0-get-str.sh. You need to provide MAPI_KEY in pd_example.py 
- generate my structure list according to FM or NM OUTCAR via 1-get-mystr.sh
- compute Ed via Ed-scanner.py
- plot via z-plot.sh

"""Curated reproducible transcription of UIUC LSAT Vol 3 Clark-Y(B) rows.
Source: CLARKYB.LFT (Cl), CLARKYB.DRG (Cd). Drag and lift sweeps
were taken separately; Cl at drag alpha is *interpolated* from the
ascending clean lift sweep and NOT a separately measured pairing.
This file defines every source value selected, so it can be audited.
"""
from pathlib import Path
import csv
import numpy as np

# Airfoil, approximately matched Reynolds groups.
# lift columns: alpha_deg, Cl; only ascending pre-stall run.
LIFT = {
  60000: [(-6.02,-.415),(-4.47,-.324),(-3.89,-.282),(-2.35,-.138),(-1.78,-.073),(-.73,.026),(.31,.125),(1.38,.258),(2.42,.352),(3.48,.498),(4.54,.632),(5.57,.708),(6.62,.804),(7.67,.901),(8.67,.984),(9.68,1.061),(10.67,1.126),(11.65,1.182),(12.62,1.210),(13.61,1.236)],
  100000: [(-6.04,-.445),(-4.39,-.296),(-3.72,-.210),(-2.56,-.080),(-1.37,.065),(-.30,.179),(.28,.263),(1.37,.385),(2.50,.513),(3.63,.627),(4.77,.741),(5.89,.849),(6.49,.905),(7.58,.994),(8.65,1.076),(9.71,1.153),(10.78,1.230),(11.85,1.265),(12.91,1.280),(13.97,1.299)],
  200000: [(-6.05,-.339),(-4.22,-.109),(-3.48,-.015),(-2.27,.105),(-1.62,.168),(-.31,.296),(.31,.369),(1.49,.501),(2.69,.614),(3.89,.725),(4.52,.785),(5.82,.900),(7.00,.999),(7.66,1.054),(8.79,1.146),(9.85,1.232),(11.00,1.288),(11.63,1.298),(12.73,1.312),(13.78,1.321)],
  300000: [(-5.93,-.241),(-4.42,-.096),(-3.43,-.009),(-2.37,.088),(-1.35,.186),(-.29,.287),(.34,.382),(1.41,.507),(2.54,.618),(3.72,.731),(4.39,.790),(5.73,.905),(7.00,1.015),(7.68,1.074),(8.86,1.177),(10.05,1.260),(10.63,1.286),(11.66,1.309),(12.71,1.331),(13.72,1.348)],
}
# drag columns: alpha_deg, Cd. Re of experiments: 60.7k, 59.5k, 100k,
# 99k, 200k, 299.8k. Second short sweep is joined to first 60k/100k.
DRAG={
  60000: [(-2.07,.0265),(-.45,.0232),(1.02,.0295),(2.54,.0333),(4.16,.0339),(5.71,.0357),(7.19,.0389),(8.74,.0390),(10.26,.0349),(11.31,.0355),(12.80,.0374)],
  100000: [(-3.07,.0249),(-1.42,.0220),(.01,.0229),(1.64,.0240),(3.17,.0225),(4.64,.0206),(6.20,.0193),(7.72,.0207),(9.29,.0255),(12.30,.0326)],
  200000: [(-4.64,.0213),(-2.92,.0156),(-1.57,.0132),(.10,.0105),(1.74,.0111),(3.11,.0115),(4.76,.0125),(6.24,.0147),(7.78,.0170),(9.30,.0192),(10.78,.0225),(12.35,.0349)],
  300000: [(-4.56,.0156),(-2.97,.0124),(-1.44,.0102),(.14,.0088),(1.71,.0090),(3.28,.0097),(4.95,.0108),(6.28,.0125),(7.87,.0148),(9.26,.0168),(10.77,.0212),(12.32,.0316)],
}

def write(out=None):
    out=Path(out or Path(__file__).with_name('clarky_uiuc_experimental.csv'))
    with out.open('w',newline='') as f:
        w=csv.writer(f)
        w.writerow(['reynolds','alpha_deg','cl_interpolated','cd_measured','cl_source','cd_source'])
        for re, dr in DRAG.items():
            alpha_lift, cl=np.array(LIFT[re]).T
            for alpha,cd in dr:
                if not (alpha_lift[0]<=alpha<=alpha_lift[-1]):
                    raise ValueError('Drag alpha outside lift coverage')
                w.writerow([re,alpha,round(float(np.interp(alpha,alpha_lift,cl)),6),cd,
                            'UIUC_LSAT_Vol3_CLARKYB.LFT','UIUC_LSAT_Vol3_CLARKYB.DRG'])
    print(out, 'rows',sum(map(len,DRAG.values())))
if __name__=='__main__':write()

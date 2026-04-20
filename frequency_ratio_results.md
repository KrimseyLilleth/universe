# Frequency-Ratio Analysis Across Scales

## Question

Given a collection of natural oscillation frequencies spanning ~28
orders of magnitude (Hubble frequency → Cs-133 hyperfine), do the
ratios between them cluster at powers of φ, 2, π, e, or any other
base more than chance predicts?

## Method

For each candidate base `b`, compute `theta_i = 2π · frac(log_b(f_i))`
and test the `theta_i` for circular uniformity using Rayleigh's R and
Kuiper's V. If the true ratios were integer powers of `b`, the
`theta_i` would concentrate at a single angle. Under the null (log-
uniform frequencies, no preferred base) they are uniform on the
circle.

Monte Carlo control: 10 000 simulated datasets of the same size, drawn
log-uniformly over the observed frequency range, to calibrate the
expected R for φ.

Frequencies compiled (N=23): Hubble, galactic-year, solar Schwabe,
solar p-mode, Earth orbit, lunar month, Earth rotation, Schumann
modes 1-4, circadian, breathing, heart rate, five EEG bands, middle C,
A4, hydrogen 21 cm, Cs-133 hyperfine. Sources are inline in
`frequency_ratio_analysis.py`.

## Result

```
base                       R   Rayleigh p   Kuiper V   Kuiper p
--------------------------------------------------------------
phi  (1.618..)        0.1443       0.6245     1.1802     0.5645
2                     0.0414       0.9622     1.0002     0.8218
e                     0.1006       0.7959     1.2202     0.5049
pi                    0.2408       0.2661     1.4102     0.2606
3                     0.3423       0.0662     1.7857     0.0400
10                    0.1830       0.4675     1.1500     0.6102
phi^2 (2.618)         0.2480       0.2453     1.4976     0.1797
sqrt(2)               0.1378       0.6509     1.3335     0.3490
```

Best raw p-value: **base 3**, Rayleigh p = 0.066.
Bonferroni-corrected for the 8 bases tested: **p = 0.53.** Not significant.

Monte Carlo: under the null, φ's R exceeds the observed value 62 % of
the time. φ beats base 2 in exactly half the null trials. The
observed φ signal is indistinguishable from noise.

Dropping the Schumann harmonics 2-4 (which are by construction
integer multiples of mode 1 and inflate any integer-ratio signal)
weakens φ further (R = 0.137) and essentially eliminates the base-2
signal (R = 0.015), while leaving the marginal base-3 hint intact.

## Conclusion

On this sample there is **no evidence** that natural oscillation
frequencies across scales cluster at powers of φ, 2, e, or π. The
only mildly anomalous base is 3, and it does not survive correction
for the number of bases tested.

Caveats, honestly stated:

1. **Sample composition matters.** Adding or removing frequencies can
   shift any of these statistics. A pre-registered list would be
   stronger.
2. **Width of each "oscillation."** Solar-cycle length varies from
   ~9 to ~14 years; EEG bands are ranges, not lines; heart rate
   depends on the subject. The point values chosen are representative
   midpoints, but the analysis treats them as exact.
3. **Harmonics inflate integer-ratio signals.** Schumann modes,
   musical octaves, and orbital resonances share mathematical
   structure by construction. The "excluding Schumann 2-4" row is
   one attempt to control for this; a thorough version would audit
   every entry.
4. **Log-uniform null may be generous.** Real frequencies cluster by
   mechanism (atomic, biological, orbital). The null could be made
   tighter with a kernel-density model of the observed distribution,
   which would lower the p-values across the board.
5. **The hypothesis could still hold within a single domain.** This
   test asks whether a single universal base governs all scales. It
   does not rule out, say, φ relations among purely biological
   oscillations or among cosmological quantities. A domain-by-domain
   version of the study is a natural follow-up.

Reproducibility: `python3 frequency_ratio_analysis.py`.

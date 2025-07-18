# Generate default report table.

    Code
      SeedMatchR::SeedMatchReport(seqs = seqs, res = res, guide.seq = guide.seq)$
        report
    Output
                      Category Full-D0 Full-D1 Full-D2 Full-D3 Full-D4 18mer-D0
      1  In silico predictions       0       1       0       0       1        0
      2  Expressed predictions       0       1       0       0       1        0
      3 Off-target predictions       0       1       0       0       1        0
      4           % off-target       0     100       0       0     100        0
        18mer-D1 18mer-D2 18mer-D3 18mer-D4 15mer-D4 15mer-D3 15mer-D2 15mer-D1
      1        0        0  3.00000 24.00000 73.00000 18.00000        0        0
      2        0        0  3.00000 24.00000 73.00000 18.00000        0        0
      3        0        0  1.00000 11.00000 21.00000  4.00000        0        0
      4        0        0 33.33333 45.83333 28.76712 22.22222        0        0
        15mer-D0 8mer 7mer-m8 7mer-A1 6mer Total           Group
      1        0    1       0       0    4 125.0 SeedMatchReport
      2        0    1       0       0    4 125.0 SeedMatchReport
      3        0    0       0       0    3  42.0 SeedMatchReport
      4        0    0       0       0   75  33.6 SeedMatchReport


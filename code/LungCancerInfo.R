LUAD.patients = list(
  Zhang = c("PA13", "PA01", "PA08", "PA17", "PA20", "PA21", "PA14", "PA09", "PA15",
            "PA11", "PA05", "PA18", "PA07", "PA04", "PA12", "PA06", "PA19", "PA16",
            "PA03"),
  Leader = c("695", "338", "370", "371", "378", "393", "403", "406", "408", "410",
             "458", "460", "464", "514", "532", "558", "564", "569", "570", "729",
             "725", "571", "572", "578", "581", "593", "596", "630", "626", "800"),
  GSE148071 = c("P2", "P5", "P8", "P9", "P12", "P13", "P16", "P20", "P21", "P24",
                "P28", "P32", "P33", "P37", "P38", "P39"),
  # GSE131907 = c("LUNG_T06","LUNG_T08","LUNG_T09","LUNG_T18","LUNG_T19","LUNG_T20",
  #               "LUNG_T30","LUNG_T31","LUNG_T34","BRONCHO_58", "EBUS_06","EBUS_28",
  #               "EBUS_49"),
  E_MTAB_6149 = c("LUAD-1", "LUAD-2"),
  GSE127465 = c("p3t", "p4t", "p5t", "p6t", "p7t"),
  GSE117570 = c("GSE117570_P1", "GSE117570_P3", "GSE117570_P4")
  # GSE153935 = c("D4T", "D5T", "D6T", "D7T", "D9T")
  # Maynard = c("TH226", "TH248", "TH179", "TH236", "TH231", "TH169", "TH158", "TH238", "AZ005", "AZ003", "AZ008")
  # Chan = c("RU1038", "RU1057", "RU1128", "RU1134", "RU1135", "RU1137", "RU1170g", "RU1262C", "RU1271")
) 
LUSC.patients = list(
  Zhang = c("PS02", "PS01", "PS08", "PS11", "PS09", "PS10", "PS07", "PS03", "PS05",
            "PS06", "PS04"),
  Leader = c("377", "522", "714", "706", "584"),
  GSE148071 = c("P1", "P3", "P4", "P6", "P7", "P10", "P14", "P15", "P17", "P18",
                "P19", "P22", "P23", "P25", "P26", "P27", "P31", "P36", "P40", "P41"),
  E_MTAB_6149 = c("LUSC-1", "LUSC-2"),
  GSE127465 = c("p1t", "p2t"),
  GSE117570 = c("GSE117570_P2")
  # GSE153935 = c("D1T", "D2T", "D3T", "D8T", "D10T", "D11T", "D12T")
)
GSE153935.LUAD = list(
  "D9T" = c("T3", "N2", "M0")
)
GSE153935.LUSC = list(
  "D3T" = c("T3", "N2", "M0"),
  "D8T" = c("T1b", "N0", "M0"),
  "D10T" = c("T1b", "N0", "M0"),
  "D11T" = c("T2a", "N0", "M0"),
  "D12T" = c("T2b", "N0", "M0")
)

# Kim = list(
#   "EBUS_06" = "IV", 
#   "EBUS_28" = "IV", 
#   "EBUS_49" = "IV", 
#   "BRONCHO_58" = "IV",
#   "LUNG_T06" = "IA", 
#   "LUNG_T08" = "IB", 
#   "LUNG_T09" = "IIA", 
#   "LUNG_T18" = "IA", 
#   "LUNG_T19" = "IA", 
#   "LUNG_T20" = "IA", 
#   "LUNG_T30" = "IA", 
#   # "LUNG_T31" = "IIIA",
#   "LUNG_T34" = "IA3"
# )

Zhang.LUAD = list(
  "PA01" = c("T1c", "N0", "M0"),
  "PA03" = c("T1b", "N0", "M0"),
  "PA04" = c("T1b", "N0", "M0"),
  "PA05" = c("T1c", "N0", "M0"),
  "PA06" = c("T1a", "N0", "M0"),
  "PA07" = c("T2a", "N0", "M0"),
  "PA08" = c("T2a", "N0", "M0"),
  "PA09" = c("T2a", "N0", "M0"),
  "PA11" = c("T2a", "N0", "M0"),
  "PA12" = c("T2a", "N0", "M0"),
  "PA13" = c("T2b", "N1", "M0"),
  "PA14" = c("T3", "N0", "M0"),
  "PA15" = c("T1b", "N2", "M0"),
  "PA16" = c("T4", "N2", "M0"),
  "PA17" = c("T3", "N2", "M0"),
  "PA18" = c("T2", "N0", "M1"),
  "PA19" = c("T3", "N2", "M1a"),
  "PA20" = c("T2a", "N2", "M1a"),
  "PA21" = c("T2b", "N0", "M1a")
)
Zhang.LUSC = list(
  "PS01" = c("T2a", "N0", "M0"),
  "PS02" = c("T2a", "N0", "M0"),
  "PS03" = c("T2a", "N0", "M0"),
  "PS04" = c("T2a", "N0", "M0"),
  "PS05" = c("T3", "N0", "M0"),
  "PS06" = c("T3", "N0", "M0"),
  "PS07" = c("T3", "N0", "M0"),
  "PS08" = c("T4", "N1", "M0"),
  "PS09" = c("T2", "N2", "M0"),
  "PS10" = c("T3", "N2", "M0"),
  "PS11" = c("T3", "N2", "M0")
)
Leader.LUAD = list(
  "338" = c("T1b", "N0", "x"),
  "370" = c("T2a", "N1", "x"),
  "378" = c("T1b", "N0", "x"),
  "393" = c("T1b", "N0", "x"),
  "403" = c("T3", "N0", "x"),
  "406" = c("T1b", "N0", "x"),
  "408" = c("T2b", "N0", "x"),
  "410" = c("T1b", "N0", "x"),
  "458" = c("T1c", "N0", "x"),
  "460" = c("T1b", "N0", "x"),
  "464" = c("T1b", "N0", "x"),
  "514" = c("T3", "N0", "x"),
  "532" = c("T1c", "N0", "x"),
  "564" = c("T1b", "N0", "x"),
  "569" = c("T1a", "N0", "x"),
  "570" = c("T2b", "N0", "x"),
  "571" = c("T1c", "N0", "x"),
  "578" = c("T1a", "N0", "x"),
  "581" = c("T2a", "N0", "x"),
  "593" = c("T1b", "N0", "x"),
  "626" = c("T1c", "N0", "x"),
  "630" = c("T1b", "N0", "x"),
  "695" = c("T2b", "N0", "x"),
  "725" = c("T1b", "N0", "x"),
  "729" = c("T2a", "N0", "x"),
  "800" = c("T1c", "N1", "x")
)
Leader.LUSC = list(
  "584" = c("T2a", "N1", "x"),
  "522" = c("T2a", "N0", "x"),
  "377" = c("T3", "N0", "x"),
  "706" = c("T2b", "N0", "x"),
  "714" = c("T2a", "N1", "NA")
)
# GSE148071 N2|N3
GSE148071.LUAD = c(
  "P2", "P5", "P8", "P9", "P12", "P13", "P16", "P20", "P21", "P24","P28", "P32", "P33", "P37", "P38", "P39"
)
GSE148071.LUSC = c(
  "P1", "P3", "P4", "P6", "P7", "P10", "P14", "P15", "P17", "P18", "P19", "P22", "P23", "P25", "P26", "P27", "P31", "P36", "P40", "P41"
)

Kim = list(
  "Y" = c("EBUS_06", "EBUS_28", "EBUS_49", "BRONCHO_58"),
  "N" = c("LUNG_T06", "LUNG_T08", "LUNG_T09", "LUNG_T18", "LUNG_T19", "LUNG_T20", "LUNG_T25", "LUNG_T30", "LUNG_T34")
)
Chan = list(
  "RU1170g" = "IIIB", 
  "RU1271" = "IV",
  "RU1128" = "IA", 
  "RU1038" = "IA"
)
Chan = list(
  "Y" = c("RU1170g", "RU1271"),
  "N" = c("RU1128", "RU1038")
)
# Maynard = list(
#   "TH158" = "IV",
#   "TH169" = "IV", 
#   "TH179" = "IV", 
#   "TH226" = "IV", 
#   "TH231" = "IV", 
#   "TH236" = "IV", 
#   "TH248" = "IIIB", 
#   "TH238" = "I",
#   "AZ005" = "IA"
# )
Maynard = list(
  "Y" = c("TH226", "TH248", "TH179", "TH236", "TH231", "TH169", "TH158"),
  "N" = c("TH238", "AZ005")
)
Lambrechts.LUAD = list(
  "LUAD-1" = c("T4", "N2", "M0"), 
  "LUAD-2" = c("T2a", "N1", "M0")
)
Lambrechts.LUSC = list(
  "LUSC-1" = c("T2b", "N0", "M0"), 
  "LUSC-2" = c("T2b", "N0", "M0")
)

# Lambrechts = list(
#   "Y" = c("LUAD-1", "LUAD-2"),
#   "N" = c("LUSC-1", "LUSC-2")
# )

Zilionis.LUAD = list(
  "p3t" = c("T4", "N0", "M1a"),
  "p4t" = c("T1b", "N0", "x"),
  "p5t" = c("T3", "N0", "M1a"),
  "p7t" = c("T3", "N1", "x")
)
Zilionis.LUSC = list(
  "p1t" = c("T2a", "N1", "x"),
  "p2t" = c("T4", "N1", "x")
)
Song.LUAD = list(
  "GSE117570_P1" = c("T1b", "N1", "x"), 
  "GSE117570_P4" = c("T1a", "N0", "x")
)
Song.LUSC = list(
  "GSE117570_P2" = c("T2a", "N0", "x")
)

# N.stage = list()
# N.stage = TNM(N.stage, Zhang)
# N.stage = TNM(N.stage, Leader)
# N.stage$Y = c(N.stage$Y, GSE148071)

LUAD.N.stage = list()
LUAD.N.stage = TNM(LUAD.N.stage, Zhang.LUAD)
LUAD.N.stage = TNM(LUAD.N.stage, Leader.LUAD)
LUAD.N.stage = TNM(LUAD.N.stage, Lambrechts.LUAD)
LUAD.N.stage = TNM(LUAD.N.stage, Zilionis.LUAD)
LUAD.N.stage = TNM(LUAD.N.stage, Song.LUAD)
# LUAD.N.stage = TNM(LUAD.N.stage, GSE153935.LUAD)
LUAD.N.stage$Y = c(LUAD.N.stage$Y, GSE148071.LUAD)

LUSC.N.stage = list()
LUSC.N.stage = TNM(LUSC.N.stage, Zhang.LUSC)
LUSC.N.stage = TNM(LUSC.N.stage, Leader.LUSC)
LUSC.N.stage = TNM(LUSC.N.stage, Lambrechts.LUSC)
LUSC.N.stage = TNM(LUSC.N.stage, Zilionis.LUSC)
LUSC.N.stage = TNM(LUSC.N.stage, Song.LUSC)
# LUSC.N.stage = TNM(LUSC.N.stage, GSE153935.LUSC)
LUSC.N.stage$Y = c(LUSC.N.stage$Y, GSE148071.LUSC)

LUAD.patients = as.character(unlist(LUAD.patients))
LUSC.patients = as.character(unlist(LUSC.patients))

# write_json(LUAD.N.stage, pretty=F, auto_unbox=TRUE, path = paste0("D:/Experiment/硕士毕业论文实验/data/LUAD.N.stage.json"))
# write_json(LUSC.N.stage, pretty=F, auto_unbox=TRUE, path = paste0("D:/Experiment/硕士毕业论文实验/data/LUSC.N.stage.json"))

# N.stage = list(
#   "Y" = c(LUAD.N.stage$Y, LUSC.N.stage$Y, Kim$Y),
#   "N" = c(LUAD.N.stage$N, LUSC.N.stage$N, Kim$N)
# )

# write_json(list("LUAD"=LUAD.patients, "LUSC"=LUSC.patients), pretty=F, auto_unbox=TRUE, path = paste0("C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/patients_subtypes.json"))
# write_json(N.stage, pretty=F, auto_unbox=TRUE, path = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/N0_vs_N1N2N3.json")

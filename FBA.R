library(sybil.tools)
library(data.table)
library(stringr)
library(sybil)
library(glpkAPI)


#Haemophilus parainfluenzae
hp1 <-readRDS("/Users/svenjabusche/Desktop/gapseq/H.parainfluenzae/H.para1.RDS")
hp2 <-readRDS("/Users/svenjabusche/Desktop/gapseq/H.parainfluenzae/H.para2.RDS")

#growthrate
sol1<- get.mod.growth(hp1)
sol2<- get.mod.growth(hp2)
sol1
sol2

#essential nutrients
dt.ess1 <- get.essential.nutrients(hp1)
dt.ess2 <- get.essential.nutrients(hp2)
dt.ess1
dt.ess2
fwrite(dt.ess1, file = "ess.Hpara_1.xls")
fwrite(dt.ess2, file = "ess.Hpara_2.xls")
View(dt.ess1)
View(dt.ess2)

#nutrient uptake
dt.nutr1 <- get.utilized.metabolites(hp1) 
dt.nutr2 <- get.utilized.metabolites(hp2) 
dt.nutr1
dt.nutr2
fwrite(dt.nutr1, file = "nutr.Hpara_1.csv")
fwrite(dt.nutr2, file = "nutr.Hpara_2.csv")
View(dt.nutr1)
View(dt.nutr2)

#produced metabolites
dt.prod1 <- get.produced.metabolites(hp1)
dt.prod2 <- get.produced.metabolites(hp2)
dt.prod1
dt.prod2
fwrite(dt.prod1, file = "prod.Hpara_1.csv")
fwrite(dt.prod2, file = "prod.Hpara_2.csv")

#exchange reactions
dt.ex1<- get.exchange.fluxes(hp1)
dt.ex2<- get.exchange.fluxes(hp1)
dt.ex1
dt.ex2
View(dt.ex1) 
View(dt.ex2) 
fwrite(dt.ex1, file = "ex.H.para1.xls")
fwrite(dt.ex2, file = "ex.H.para2.xls")

#Streptococcus vestibularis
sv1 <-readRDS("/Users/svenjabusche/Desktop/gapseq/S.vestibularis/S.vestibularis_1.RDS")
sv2 <-readRDS("/Users/svenjabusche/Desktop/gapseq/S.vestibularis/S.vestibularis_2.RDS")

#growthrate
sol3<- get.mod.growth(sv1)
sol4<- get.mod.growth(sv2)
sol3
sol4

#essential nutrients
dt.ess3 <- get.essential.nutrients(sv1)
dt.ess4 <- get.essential.nutrients(sv2)
dt.ess3
dt.ess4
View(dt.ess3)
View(dt.ess4)

#nutrient uptake
dt.nutr3 <- get.utilized.metabolites(sv1) 
dt.nutr4 <- get.utilized.metabolites(sv2) 
dt.nutr3
dt.nutr4
View(dt.nutr3)
View(dt.nutr4)
fwrite(dt.nutr3, file = "nutr.Svest_1.csv")
fwrite(dt.nutr4, file = "nutr.Svest_2.csv")

#produced metabolites
dt.prod3 <- get.produced.metabolites(sv1)
dt.prod4 <- get.produced.metabolites(sv2)
dt.prod3
dt.prod4
fwrite(dt.prod3, file = "prod.Svest_1.csv")
fwrite(dt.prod4, file = "prod.Svest_2.csv")

#exchange reactions
dt.ex3<- get.exchange.fluxes(sv1)
dt.ex4<- get.exchange.fluxes(sv2)
View(dt.ex3) 
View(dt.ex4) 
fwrite(dt.ex3, file = "ex.Svest1.xls")
fwrite(dt.ex4, file = "ex.Svest2.xls")


#Rothia mucilaginosa
rm1 <-readRDS("/Users/svenjabusche/Desktop/gapseq/Rothia/R.mucila_1.RDS")
rm2 <-readRDS("/Users/svenjabusche/Desktop/gapseq/Rothia/R.mucila_2.RDS")

#growthrate
sol5<- get.mod.growth(rm1)
sol6<- get.mod.growth(rm2)
sol5
sol6

#essential nutrients
dt.ess5 <- get.essential.nutrients(rm1)
dt.ess6 <- get.essential.nutrients(rm2)
dt.ess5
dt.ess6
fwrite(dt.ess5, file = "ess.Rmuci_1.xls")
fwrite(dt.ess6, file = "ess.Rmuci_2.xls")
View(dt.ess5)
View(dt.ess6)

#nutrient uptake
dt.nutr5 <- get.utilized.metabolites(rm1) 
dt.nutr6 <- get.utilized.metabolites(rm2) 
dt.nutr5
dt.nutr6
View(dt.nutr5)
View(dt.nutr6)
fwrite(dt.nutr5, file = "nutr.Rmuci_1.csv")
fwrite(dt.nutr6, file = "nutr.Rmuci_2.csv")

#produced metabolites
dt.prod5 <- get.produced.metabolites(rm1)
dt.prod6 <- get.produced.metabolites(rm2)
dt.prod5
dt.prod6
fwrite(dt.prod5, file = "prod.Rmuci_1.csv")
fwrite(dt.prod6, file = "prod.Rmuci_2.csv")

#exchange reactions
dt.ex5<- get.exchange.fluxes(rm1)
dt.ex6<- get.exchange.fluxes(rm2)
View(dt.ex5) 
View(dt.ex6) 
fwrite(dt.ex5, file = "ex.Rmuci_1.xls")
fwrite(dt.ex6, file = "ex.Rmuci_2.xls")


#Prevotella melaninogenica
pm1 <-readRDS("/Users/svenjabusche/Desktop/gapseq/P.mela/P.mela1.RDS")
pm2 <-readRDS("/Users/svenjabusche/Desktop/gapseq/P.mela/P.mela2.RDS")

#growthrate
sol7<- get.mod.growth(pm1)
sol8<- get.mod.growth(pm2)
sol7
sol8

#essential nutrients
dt.ess7 <- get.essential.nutrients(pm1)
dt.ess8 <- get.essential.nutrients(pm2)
dt.ess7
dt.ess8
fwrite(dt.ess7, file = "ess.P.mela1.csv")
fwrite(dt.ess8, file = "ess.P.mela2.csv")
View(dt.ess7)
View(dt.ess8)

#nutrient uptake
dt.nutr7 <- get.utilized.metabolites(pm1) 
dt.nutr8 <- get.utilized.metabolites(pm2) 
dt.nutr7
dt.nutr8
fwrite(dt.nutr7, file = "nutr.P.mela1.csv")
fwrite(dt.nutr8, file = "nutr.P.mela2.csv")

#produced metabolites
dt.prod7 <- get.produced.metabolites(pm1)
dt.prod8 <- get.produced.metabolites(pm2)
dt.prod7
dt.prod8
fwrite(dt.prod7, file = "prod.P.mela1.csv")
fwrite(dt.prod8, file = "prod.P.mela2.csv")

#exchange reactions
dt.ex7 <- get.exchange.fluxes(pm1)
dt.ex8 <- get.exchange.fluxes(pm2)
View(dt.ex7) 
View(dt.ex8) 
fwrite(dt.ex8, file = "ex.P.mela1.csv")
fwrite(dt.e82, file = "ex.P.mela2.csv")


#Veilonella parvula
vp1 <-readRDS("/Users/svenjabusche/Desktop/gapseq/V.parv/V.parv_1.RDS")
vp2 <-readRDS("/Users/svenjabusche/Desktop/gapseq/V.parv/V.parv_2.RDS")

#growthrate
sol9<- get.mod.growth(vp1)
sol10<- get.mod.growth(vp2)
sol9
sol10

#essential nutrients
dt.ess9 <- get.essential.nutrients(vp1)
dt.ess10 <- get.essential.nutrients(vp2)
dt.ess9
dt.ess10
fwrite(dt.ess9, file = "ess.V.parv1.csv")
fwrite(dt.ess10, file = "ess.V.parv2.csv")
View(dt.ess9)
View(dt.ess10)

#nutrient uptake
dt.nutr9 <- get.utilized.metabolites(vp1) 
dt.nutr10 <- get.utilized.metabolites(vp2) 
dt.nutr9
dt.nutr10
fwrite(dt.nutr9, file = "nutr.V.parv1.csv")
fwrite(dt.nutr10, file = "nutr.V.parv2.csv")

#produced metabolites
dt.prod9 <- get.produced.metabolites(vp1)
dt.prod10 <- get.produced.metabolites(vp2)
dt.prod9
dt.prod10
fwrite(dt.prod9, file = "prod.V.parv1.csv")
fwrite(dt.prod10, file = "prod.V.parv2.csv")

#exchange reactions
dt.ex9 <- get.exchange.fluxes(vp1)
dt.ex10 <- get.exchange.fluxes(vp2)
View(dt.ex9) 
View(dt.ex10) 
fwrite(dt.ex9, file = "ex.V.parv1.csv")
fwrite(dt.ex10, file = "ex.V.parv2.csv")


#Fusobacterium periodonticum
fp1 <-readRDS("/Users/svenjabusche/Desktop/gapseq/F.peri/F.peri1.RDS")
fp2 <-readRDS("/Users/svenjabusche/Desktop/gapseq/F.peri/F.peri2.RDS")

#growthrate
sol11<- get.mod.growth(fp1)
sol12<- get.mod.growth(fp2)
sol11
sol12

#essential nutrients
dt.ess11 <- get.essential.nutrients(fp1)
dt.ess12 <- get.essential.nutrients(fp2)
dt.ess11
dt.ess12
fwrite(dt.ess11, file = "ess.F.peri1.csv")
fwrite(dt.ess12, file = "ess.F.peri2.csv")
View(dt.ess11)
View(dt.ess12)

#nutrient uptake
dt.nutr11 <- get.utilized.metabolites(fp1) 
dt.nutr12 <- get.utilized.metabolites(fp2) 
dt.nutr11
dt.nutr12
fwrite(dt.nutr11, file = "nutr.F.peri1.csv")
fwrite(dt.nutr12, file = "nutr.F.peri2.csv")

#produced metabolites
dt.prod11 <- get.produced.metabolites(fp1)
dt.prod12 <- get.produced.metabolites(fp2)
dt.prod11
dt.prod12
fwrite(dt.prod11, file = "prod.F.peri1.csv")
fwrite(dt.prod12, file = "prod.F.peri2.csv")

#exchange reactions
dt.ex11 <- get.exchange.fluxes(fp1)
dt.ex12 <- get.exchange.fluxes(fp2)
View(dt.ex11) 
View(dt.ex12) 
fwrite(dt.ex11, file = "ex.F.peri1.csv")
fwrite(dt.ex12, file = "ex.F.peri2.csv")

                       
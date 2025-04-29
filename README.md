# bioinformatica2025
# Consultar la pagina de markdown para saber los codigos

https://www.markdownguide.org/basic-syntax/

## para hacer la anova en r consultar este link: https://statsandr.com/blog/anova-in-r/

##### imagen de la anova
![unnamed-chunk-12-1](https://github.com/user-attachments/assets/205a303c-943d-45de-8f36-9aca2a05b6e9)

# script genomica de poblaciones

2025.03.31
##########################################Microsatélites#####################################

setwd("C:/Users/Alumno/Desktop")
 library(adegenet)
obj <- read.genetix("datapopR.gtx")
catpop <- genind2genpop(obj)
 locNames(catpop)
 class(obj)
 summary(obj)
library(pegas)
 round(pegas::hw.test(obj,B=100),digits = 3)
HWE.test <- data.frame(sapply(seppop(obj), 
                               function(ls) pegas::hw.test(ls, B=0)[,3]))
 HWE.test.chisq <- t(data.matrix(HWE.test))
 {cat("Chi-squared test (p-values):", "\n")
     round(HWE.test.chisq,3)}
 HWE.test <- data.frame(sapply(seppop(obj), 
                               function(ls) pegas::hw.test(ls, B=100)[,4]))
 HWE.test.MC <- t(data.matrix(HWE.test))
 {cat("MC permuation test (p-values):", "\n")
     round(HWE.test.MC,3)}
 alpha=0.05
 Prop.loci.out.of.HWE <- data.frame(Chisq=apply(HWE.test.chisq<alpha, 2, mean), 
                                    MC=apply(HWE.test.MC<alpha, 2, mean))
 Prop.loci.out.of.HWE     
 Prop.pops.out.of.HWE <- data.frame(Chisq=apply(HWE.test.chisq<alpha, 1, mean), 
                                    MC=apply(HWE.test.MC<alpha, 1, mean))
 Prop.pops.out.of.HWE  
Chisq.fdr <- matrix(p.adjust(HWE.test.chisq,method="fdr"), 
                    nrow=nrow(HWE.test.chisq))
MC.fdr <- matrix(p.adjust(HWE.test.MC, method="fdr"), 
                    nrow=nrow(HWE.test.MC))

Prop.pops.out.of.HWE <- data.frame(Chisq=apply(HWE.test.chisq<alpha, 1, mean), 
           MC=apply(HWE.test.MC<alpha, 1, mean),
           Chisq.fdr=apply(Chisq.fdr<alpha, 1, mean),
           MC.fdr=apply(MC.fdr<alpha, 1, mean))
Prop.pops.out.of.HWE 
 
#install.packages("poppr")
 library(poppr)
poppr::ia(obj, sample=199) 


LD.pair <- poppr::pair.ia(obj)

#install.packages("PopGenReport")
library(PopGenReport)

Null.alleles <- PopGenReport::null.all(obj)
Null.alleles

Null.alleles$homozygotes$probability.obs.    




###########################################SNPS###################
##https://datadryad.org/dataset/doi:10.5061/dryad.0rxwdbrwf
#Rapid and strong population genetic differentiation and genomic signatures of climatic adaptation in an invasive mealybug





# Cargar librerías necesarias
library(vcfR)
library(adegenet)
library(pegas)
library(poppr)
library(PopGenReport)

# Leer archivo VCF

vcf <- read.vcfR("all.SNP.ps.11pop.19634.vcf")

# Convertir a genind (funciona bien para SNPs con alelos binarios)
obj <- vcfR2genind(vcf, sep = "/", ncode = 1)

# Asignar poblaciones (primeros 4 caracteres de nombre del individuo)
pop(obj) <- substr(indNames(obj), 1, 4)

# Verificar agrupamiento
table(pop(obj))

######continuar script normal como los micros#######
catpop <- genind2genpop(obj)
 locNames(catpop)
 class(obj)
 summary(obj)
library(pegas)
 round(pegas::hw.test(obj,B=100),digits = 3)
HWE.test <- data.frame(sapply(seppop(obj), 
                               function(ls) pegas::hw.test(ls, B=0)[,3]))
 HWE.test.chisq <- t(data.matrix(HWE.test))
 {cat("Chi-squared test (p-values):", "\n")
     round(HWE.test.chisq,3)}
 HWE.test <- data.frame(sapply(seppop(obj), 
                               function(ls) pegas::hw.test(ls, B=100)[,4]))
 HWE.test.MC <- t(data.matrix(HWE.test))
 {cat("MC permuation test (p-values):", "\n")
     round(HWE.test.MC,3)}
 alpha=0.05
 Prop.loci.out.of.HWE <- data.frame(Chisq=apply(HWE.test.chisq<alpha, 2, mean), 
                                    MC=apply(HWE.test.MC<alpha, 2, mean))
 Prop.loci.out.of.HWE     
 Prop.pops.out.of.HWE <- data.frame(Chisq=apply(HWE.test.chisq<alpha, 1, mean), 
                                    MC=apply(HWE.test.MC<alpha, 1, mean))
 Prop.pops.out.of.HWE  
Chisq.fdr <- matrix(p.adjust(HWE.test.chisq,method="fdr"), 
                    nrow=nrow(HWE.test.chisq))
MC.fdr <- matrix(p.adjust(HWE.test.MC, method="fdr"), 
                    nrow=nrow(HWE.test.MC))

Prop.pops.out.of.HWE <- data.frame(Chisq=apply(HWE.test.chisq<alpha, 1, mean), 
           MC=apply(HWE.test.MC<alpha, 1, mean),
           Chisq.fdr=apply(Chisq.fdr<alpha, 1, mean),
           MC.fdr=apply(MC.fdr<alpha, 1, mean))
Prop.pops.out.of.HWE 
 
#install.packages("poppr")
 library(poppr)
poppr::ia(obj, sample=199) 


LD.pair <- poppr::pair.ia(obj)

#install.packages("PopGenReport")
library(PopGenReport)

Null.alleles <- PopGenReport::null.all(obj)
Null.alleles

Null.alleles$homozygotes$probability.obs.   




#######################SNPS_ REDUCIDO_CLASE####################################
## --------------------------------------------
# 01_preparar_datos_SNPs.R
# Autor: Armando Sunny
# Objetivo: Cargar archivo VCF, filtrar poblaciones y submuestrear loci polimórficos
# --------------------------------------------

# -------------------------
# 1. Cargar librerías
# -------------------------
library(vcfR)
library(adegenet)
library(pegas)
library(poppr)
library(PopGenReport)

# -------------------------
# 2. Leer archivo VCF
# -------------------------
vcf <- read.vcfR("all.SNP.ps.11pop.19634.vcf")

# Convertir a objeto genind
obj <- vcfR2genind(vcf, sep = "/", ncode = 1)

# Asignar poblaciones (primeros 4 caracteres del nombre del individuo)
pop(obj) <- substr(indNames(obj), 1, 4)

# -------------------------
# 3. Filtrar poblaciones grandes (>= 18 individuos)
# -------------------------
pop_sizes <- table(pop(obj))
pops_to_keep <- names(pop_sizes[pop_sizes >= 18])

obj_filt <- obj[pop(obj) %in% pops_to_keep, ]
cat("Poblaciones incluidas:\n")
print(table(pop(obj_filt)))

# -------------------------
# 4. Submuestreo simple de 300 loci (para clase)
# -------------------------
set.seed(123)
n_snp <- 300
obj_small <- obj_filt[, sample(1:nLoc(obj_filt), n_snp)]
locNames(obj_small) <- paste0("Locus_", seq_along(locNames(obj_small)))

# Guardar para usar en clase
save(obj_small, file = "SNPs_para_clase_submuestreados.RData")
cat("Objeto 'obj_small' guardado: SNPs_para_clase_submuestreados.RData\n")

# -------------------------
# 5. CREAR VERSIÓN SOLO CON LOCI POLIMÓRFICOS
# -------------------------
# Identificar todos los loci polimórficos en el objeto original (antes del filtrado de poblaciones)
all_poly <- which(obj@loc.n.all > 1)

# Verificar cuántos loci hay disponibles
cat("Total de loci polimórficos en obj original:", length(all_poly), "\n")

# Submuestrear hasta 300 loci polimórficos (o menos si no hay suficientes)
set.seed(123)
selected_loci <- sample(all_poly, min(300, length(all_poly)))

# Crear objeto con loci polimórficos
obj_polyonly <- obj[, selected_loci]
locNames(obj_polyonly) <- paste0("Locus_", seq_along(locNames(obj_polyonly)))

# Revisar
cat("Loci en obj_polyonly:", length(locNames(obj_polyonly)), "\n")

# Guardar para análisis
save(obj_polyonly, file = "SNPs_para_clase_POLIMORFICOS.RData")
cat("Objeto 'obj_polyonly' guardado: SNPs_para_clase_POLIMORFICOS.RData\n")












######continuar script normal como los micros#######
catpop <- genind2genpop(obj_polyonly)
 locNames(catpop)
 class(obj_polyonly)
 summary(obj_polyonly)
library(pegas)



##########################Hardy-Weinberg ##########################

 round(pegas::hw.test(obj_polyonly,B=100),digits = 3)
HWE.test <- data.frame(sapply(seppop(obj_polyonly), 
                               function(ls) pegas::hw.test(ls, B=0)[,3]))
 HWE.test.chisq <- t(data.matrix(HWE.test))
 {cat("Chi-squared test (p-values):", "\n")
     round(HWE.test.chisq,3)}
 HWE.test <- data.frame(sapply(seppop(obj_polyonly), 
                               function(ls) pegas::hw.test(ls, B=100)[,4]))
 HWE.test.MC <- t(data.matrix(HWE.test))
 {cat("MC permuation test (p-values):", "\n")
     round(HWE.test.MC,3)}
 alpha=0.05
 Prop.loci.out.of.HWE <- data.frame(Chisq=apply(HWE.test.chisq<alpha, 2, mean), 
                                    MC=apply(HWE.test.MC<alpha, 2, mean))
 Prop.loci.out.of.HWE     
 Prop.pops.out.of.HWE <- data.frame(Chisq=apply(HWE.test.chisq<alpha, 1, mean), 
                                    MC=apply(HWE.test.MC<alpha, 1, mean))
 Prop.pops.out.of.HWE  
Chisq.fdr <- matrix(p.adjust(HWE.test.chisq,method="fdr"), 
                    nrow=nrow(HWE.test.chisq))
MC.fdr <- matrix(p.adjust(HWE.test.MC, method="fdr"), 
                    nrow=nrow(HWE.test.MC))

Prop.pops.out.of.HWE <- data.frame(Chisq=apply(HWE.test.chisq<alpha, 1, mean), 
           MC=apply(HWE.test.MC<alpha, 1, mean),
           Chisq.fdr=apply(Chisq.fdr<alpha, 1, mean),
           MC.fdr=apply(MC.fdr<alpha, 1, mean))
Prop.pops.out.of.HWE 
 




############ Hardy-Weinberg por población con manejo de errores#################
HWE.test.chisq <- sapply(seppop(obj_polyonly), function(ls) {
  tryCatch({
    pegas::hw.test(ls, B = 0)[,3]  # valor p del test Chi-cuadrado
  }, error = function(e) {
    message("Error en hw.test: ", e$message)
    rep(NA, nLoc(ls))
  })
})

# Transponer y convertir a matriz
HWE.test.chisq <- t(data.matrix(HWE.test.chisq))

# Ver resultados redondeados
round(HWE.test.chisq, 3)


# Proporción de loci fuera de HWE por población
alpha <- 0.05
Prop.pops.out.of.HWE <- data.frame(
  Chisq = apply(HWE.test.chisq < alpha, 1, mean, na.rm = TRUE)
)

Prop.pops.out.of.HWE


###########ALELOS_NULOS############################



#install.packages("poppr")
 library(poppr)
poppr::ia(obj_polyonly, sample=199) 


LD.pair <- poppr::pair.ia(obj_polyonly)

#install.packages("PopGenReport")
library(PopGenReport)

Null.alleles <- PopGenReport::null.all(obj_polyonly)
Null.alleles

Null.alleles$homozygotes$probability.obs.   


{cat(" summary1 (Chakraborty et al. 1994):", "\n")
round(Null.alleles$null.allele.freq$summary1,2)} 

{cat("summary2 (Brookfield et al. 1996):", "\n")
round(Null.alleles$null.allele.freq$summary2,2)}   


########diversidad Genetica###############


Sum <- adegenet::summary(obj_polyonly)
names(Sum)

Sum



par(mar=c(5.5, 4.5,1,1))
barplot(Sum$pop.n.all, las=3, 
       xlab = "", ylab = "Number of alleles")



plot(Sum$n.by.pop, Sum$pop.n.all, 
       xlab = "Sample size", ylab = "Number of alleles")
abline(lm(Sum$pop.n.all ~ Sum$n.by.pop), col = "red")






#install.packages("hierfstat")  # si no está instalada
library(hierfstat)





# Convertir a formato hierfstat
obj_hier <- genind2hierfstat(obj_polyonly)

# Asegurar que la columna 'pop' está presente
obj_hier$pop <- as.numeric(pop(obj_polyonly))

# Calcular rarefied allelic richness
Richness <- allelic.richness(obj_hier)

# Calcular promedio de riqueza por población
mean.rich <- rowMeans(Richness$Ar, na.rm = TRUE)

# Verificar valores calculados
cat("\nRiqueza alélica por población (media rarificada):\n")
print(round(mean.rich, 3))

# Gráfico de barras
par(mar = c(6, 5, 2, 1))
barplot(mean.rich,
        las = 3,
        ylab = "Rarefied allelic richness (Ar)",
        col = "lightblue",
        main = "Riqueza alélica rarificada por población")

# Opcional: guardar como imagen
# png("riqueza_alelica.png", width = 800, height = 600)
# barplot(mean.rich, las = 3, ylab = "Rarefied allelic richness (Ar)", col = "skyblue")
# dev.off()


  par(mar=c(3, 4.5,1,1))
  barplot(Sum$Hexp, ylim=c(0,1), ylab="Expected heterozygosity")
  barplot(Sum$Hobs, ylim=c(0,1), ylab="Observed heterozygosity")






  Hobs <- t(sapply(seppop(obj_small), function(ls) summary(ls)$Hobs))
  Hexp <- t(sapply(seppop(obj_small), function(ls) summary(ls)$Hexp))
  {cat("Expected heterozygosity (Hexp):", "\n")
  round(Hexp, 2)}


  {cat("\n", "Observed heterozygosity (Hobs):", "\n")
  round(Hobs, 2)}





  par(mar=c(5.5, 4.5, 1, 1))
  Hobs.pop <- apply(Hobs, MARGIN = 1, FUN = mean)
  Hexp.pop <- apply(Hexp, 1, mean) 
  barplot(Hexp.pop, ylim=c(0,1), las=3, ylab="Expected heterozygosity")
  barplot(Hobs.pop, ylim=c(0,1), las=3, ylab="Observed heterozygosity")




##Aggregate genetic data at population level (allele frequencies)


Frogs.genpop <- adegenet::genind2genpop(obj_polyonly)


Freq <- adegenet::makefreq(Frogs.genpop)


round(Freq[1:6,1:10], 2)


apply(Freq, MARGIN = 1, FUN = sum)    # Just checking



########ESTRUCTURA_GENETICA#######

# Revisar poblaciones
table(pop(obj_polyonly))

# 03_amova.R - Análisis AMOVA con jerarquía Cluster/POP

# ----------------------
# 1. Cargar librerías
# ----------------------
library(adegenet)
library(poppr)
library(pegas)
library(ade4)

# ----------------------
# 2. Cargar objeto SNPs submuestreados
# ----------------------
load("SNPs_para_clase_submuestreados.RData")  # Contiene `obj_polyonly`

# ----------------------
# 3. Asignar clústeres por población
# ----------------------
cluster_assign <- c(
  "AHHF" = "Norte", "FJPT" = "Norte",
  "GDGZ" = "Centro", "GXNN" = "Centro", "HNCS" = "Centro", 
  "HNYY" = "Centro", "JSWX" = "Centro",
  "JXGZ" = "Sur", "JXNC" = "Sur", "XJTL" = "Sur", "ZJHZ" = "Sur"
)

# ----------------------
# 4. Crear data.frame con strata
# ----------------------
df <- genind2df(obj_polyonly, sep = "/")
df$POP <- pop(obj_polyonly)
df$Cluster <- cluster_assign[as.character(df$POP)]

# ----------------------
# 5. Reconstruir objeto genind con jerarquía
# ----------------------
obj_poly_hier <- df2genind(
  X = df[, !(names(df) %in% c("POP", "Cluster"))],
  sep = "/", ncode = 1,
  ploidy = 2,
  pop = df$POP,
  strata = df[, c("Cluster", "POP")],
  type = "codom"
)



# ----------------------
#6. Checar que esten las poblaciones 
# ----------------------


table(df$POP)  # Poblaciones antes de reconstruir el objeto

table(pop(obj_poly_hier))  # Poblaciones después de reconstrucción

unique(df$POP)[!(unique(df$POP) %in% names(cluster_assign))]





# ----------------------
# 7. Ejecutar AMOVA
# ----------------------
amova_result <- poppr.amova(obj_poly_hier, hier = ~Cluster/POP, method = "ade4")
print(amova_result)

# ----------------------
# 8. Permutación para significancia
# ----------------------
set.seed(999)
amova_test <- ade4::randtest(amova_result, nrepet = 999)

# ----------------------
# 9. Guardar resultados
# ----------------------
write.csv(as.data.frame(amova_result$results), "AMOVA_Resultados.csv")
pdf("AMOVA_permutation_plot.pdf")
plot(amova_test)
dev.off()

# ----------------------
# 10. Mensaje final
# ----------------------
cat("\n\u2714\ufe0f AMOVA y test de permutación completados. Resultados exportados.\n")





########PCA###################################
####

#install.packages("dartR.base")
#BiocManager::install("SNPRelate")
#BiocManager::install("snpStats")

library(adegenet)
library(methods)
library(dartR.base)
library(SNPRelate)
library(snpStats)


# Convertir genind a genlight
genlight <- gi2gl(obj_polyonly)




# Asignar poblaciones
pop(genlight) <- pop(obj_polyonly)

# Eliminar loci con demasiados datos faltantes
toRemove <- is.na(glMean(genlight, alleleAsUnit = TRUE))
genlight <- genlight[, !toRemove]

# Detectar estructura genética (agrupamiento)
cluster <- find.clusters(genlight,
                         max.n.clust = 11,
                         stat = "BIC",
                         n.iter = 1e6,
                         n.start = 1000,
                         scale = FALSE,
                         parallel = TRUE)
# Ejecutar DAPC con el número de PCs y discriminantes que elijas
dapc_result <- dapc(genlight, pop = cluster$grp, n.pca = 50, n.da = 2)


# 8. Visualización
# Dispersión en el espacio discriminante
scatter(dapc_result, scree.da = TRUE, scree.pca = TRUE)

# Composición de grupos
compoplot(dapc_result, show.lab = TRUE, col = funky(3), legend = TRUE)

# 9. Exportar resultados
write.csv(as.data.frame(dapc_result$posterior), file = "DAPC_grupos_probabilidades.csv")
write.csv(data.frame(Individuo = rownames(dapc_result$posterior), Grupo = cluster$grp),
          file = "DAPC_asignacion_individuos.csv")












#######Calcular F-statistics (FST, FIS, FIT) con hierfstat#####################



library(hierfstat)

# Convertir a formato hierfstat
obj_hier <- genind2hierfstat(obj_polyonly)

# Calcular F-statistics globales y por locus
fstats <- hierfstat::basic.stats(obj_hier)

# Valores globales (por población)
fstats$overall



#D de Jost (con mmod)


jost_d <- genet.dist(obj_polyonly, method = "D")

jost_d

#Exportar resultados a CSV (opcional)

# Fst, Fis por locus
write.csv(fstats$perloc, "FST_por_locus.csv")

# Jost D
write.csv(as.matrix(jost_d), "D_jost_matrix.csv")




# -----------------------------------------------
# SCRIPT AJUSTADO PARA DATOS DE PARENTESCO CON SNPRelate
# -----------------------------------------------

# 1. Instalar y cargar librería necesaria
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SNPRelate")

library(SNPRelate)

# 2. Definir nombres de archivo
vcf.fn <- "all.SNP.ps.11pop.19634.vcf"
gds.fn <- "genotipo.gds"

# 3. Convertir archivo VCF a GDS
snpgdsVCF2GDS(vcf.fn, gds.fn, method = "biallelic.only")

# 4. Abrir archivo GDS
genofile <- snpgdsOpen(gds.fn)

# --------------------------
# OPCIÓN A: MATRIZ DE IBS (no requiere información de cromosomas)
# --------------------------

# 5. Cargar objeto con los individuos seleccionados
load("SNPs_para_clase_POLIMORFICOS.RData")  # carga obj_polyonly

# 6. Obtener IDs de individuos seleccionados
ind_keep <- indNames(obj_polyonly)

# 7. Ejecutar IBS solo para esos individuos
ibs <- snpgdsIBS(genofile, sample.id = ind_keep, autosome.only = FALSE)

# Extraer matriz
ibs_mat <- ibs$ibs

# Ver dimensiones y visualizar
dim(ibs_mat)
heatmap(ibs_mat, symm = TRUE, Rowv = NA, Colv = NA,
        main = "Matriz de identidad por estado (IBS)")

# Exportar a CSV
write.csv(ibs_mat, file = "IBS_matrix_filtrada.csv")

# --------------------------
# EXTRA: Si alguna vez necesitas IBD (requiere cromosomas válidos en VCF)
# --------------------------
# ibd <- snpgdsIBDMLE(genofile)  # <- Esta línea falla si no hay cromosomas válidos
# heatmap(as.matrix(ibd$kinship), Rowv=NA, Colv=NA)

# --------------------------
# Cierre
# --------------------------
snpgdsClose(genofile)






#######MIGRACION_RECIENTE#############

assignplot(dapc_result)



#####parentesco y estructura familiar dendograma######



hc <- hclust(as.dist(1 - ibs_mat), method = "average")
plot(hc, main = "Clustering por IBS")






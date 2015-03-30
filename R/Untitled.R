library(Momocs)
load("~/Desktop/allF.rda")

allP <- PCA(allF)
MANOVA_PW(filter(allP, !Island %in% c("New_Zealand", "Pitcairn")), "Island")


New_Zealand_scores <- filter(allP, Island=="New_Zealand")$x[, 1]
New_Britain_scores <- filter(allP, Island=="New_Britain")$x[, 1]
Chatham_scores <- filter(allP, Island=="Chatham")$x[, 1]
Pitcairn_scores <- filter(allP, Island=="Pitcairn")$x[, 1]
Rapa_Nui_scores <- filter(allP, Island=="Rapa_Nui")$x[, 1]

wilcox.test(New_Zealand_scores, New_Britain_scores)
wilcox.test(New_Zealand_scores, Chatham_scores)
wilcox.test(New_Zealand_scores, Pitcairn_scores)
wilcox.test(New_Zealand_scores, Rapa_Nui_scores)

wilcox.test(Chatham_scores, New_Britain_scores)
wilcox.test(Chatham_scores, New_Zealand_scores)
wilcox.test(Chatham_scores, Pitcairn_scores)
wilcox.test(Chatham_scores, Rapa_Nui_scores)


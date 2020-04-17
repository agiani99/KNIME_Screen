library(tidyverse)
library(readxl)
library(caret)
library(xgboost)
#library(classifierplots)
library(knitr)
library(broom)

load("W:/projects/ESP Projects Active/Corona_Broad_2020/GeneratedData/DataForDistribution/Figures_&_Tables/20200405_PS_all_valid_v3.RData")

ALL_valid_v3 <- read_excel("20200331_PS_all_valid_v3.xlsx", 
                           col_types = c("text", "skip", "text", 
                                         "text", "text", "text", "text", "text", 
                                         "text", "text", "skip", "skip", "numeric", 
                                         "numeric", "text", "skip", "skip", 
                                         "skip", "numeric", "numeric", "numeric", 
                                         "numeric", "numeric", "numeric", 
                                         "numeric", "numeric", "numeric", 
                                         "numeric", "numeric", "numeric", 
                                         "numeric", "numeric", "numeric", 
                                         "numeric", "numeric", "numeric", 
                                         "numeric", "numeric", "numeric", 
                                         "numeric", "numeric", "numeric", 
                                         "numeric", "skip", "text", "skip", 
                                         "text", "numeric", "numeric", "numeric", 
                                         "numeric", "numeric", "numeric", 
                                         "numeric", "numeric", "numeric", 
                                         "numeric", "numeric", "numeric", 
                                         "numeric", "numeric", "numeric", 
                                         "numeric", "numeric", "numeric", 
                                         "text"))


names(ALL_valid_v3)

ALL_valid_v3$abs_ID <- 1:nrow(ALL_valid_v3)

ALL_valid_num <- ALL_valid_v3 %>% select_if(.,list(is.numeric)) %>% select(abs_ID, everything())

names(ALL_valid_num) <- gsub(" ", "", names(ALL_valid_num))
names(ALL_valid_v3) <- gsub(" ", "", names(ALL_valid_v3))

tidyselect::vars_select(names(ALL_valid_num), matches('?[\\)]'))

#View(ALL_valid_num)

ALL_valid_num_df  <- ALL_valid_v3 %>% filter(Column != 12) %>% 
  select(c("abs_ID", tidyselect::vars_select(names(ALL_valid_num), matches('?[\\)]')))) %>% #names()
  select(-c("%Confluence(DPC)", "NumberOfCells(DPC)", "%Inhibition(NucleiStaining)",
            "NumberofAnalyzedFields(NucleiStaining)" ))

ALL_valid_num_df$Active75 <- if_else(ALL_valid_num_df$`%Inhibition(DPC)` < 75, "inactive75", "active75")

ALL_valid_num_df <- ALL_valid_num_df %>% select(abs_ID, Active75, everything()) %>% select(-`%Inhibition(DPC)`)

ALL_valid_num_df$Active75 <- if_else(ALL_valid_num_df$Active75 == "active75", 1,0)

############################################################
#                                                          #
#                      xgboost model                       #
#                                                          #
############################################################

set.seed(2401)

trainIndex <- as.vector(createDataPartition(ALL_valid_num_df$Active75, p = .80, 
                                  list = FALSE, 
                                  times = 1))

MeTrain <- ALL_valid_num_df[trainIndex,]
MeTest  <- ALL_valid_num_df[-trainIndex ,]
sum(is.na(MeTrain[,-c(1:2)]))

labeltrain <- MeTrain$Active75
labeltest <- MeTest$Active75
table(labeltest)
table(labeltrain)
#tail(names(MeTrain))
#head(names(MeTrain))
#str(MeTest)

dtrain <- xgb.DMatrix(data = as.matrix(MeTrain[,-c(1,2)]), label = labeltrain) 
dtest <- xgb.DMatrix(data = as.matrix(MeTest[,-c(1,2)]),label = labeltest)

xgb.fit <- xgb.train(data = dtrain, label = labeltrain, booster = "gbtree",
                     objective = "binary:logistic",
                     colsample_bytree = 0.5, gamma = 0,
                     learning_rate = 0.05, max_depth = 5,
                     min_child_weight = 1,
                     #reg_alpha = 0.9, reg_lambda = 0.4,
                     subsample = 0.9, print_every_n = 1,
                     silent = 1, nrounds = 100, watchlist=list(train=dtrain, test=dtest))

xgb.fit

cv.xgb <- xgb.cv(data = dtrain, nrounds = 500, nfold = 10, label = labeltrain,
                 objective = "binary:logistic", #booster = "gblinear",
                 colsample_bytree = 0.5, gamma = 0, booster = "gbtree",
                 learning_rate = 0.007, max_depth = 6,
                 min_child_weight = 1, #n_estimators = 7300,
                 reg_alpha = 0.7, reg_lambda = .2, print_every_n = 1,
                 subsample = 0.9)
cv.xgb

xgb.best <- xgboost(data = dtrain, label = labeltrain,
                    booster = "gbtree", objective = "binary:logistic",
                    colsample_bytree = 0.5, gamma = 0,
                    learning_rate = 0.01, max_depth = 5,
                    min_child_weight = 1, #n_estimators = 7300,
                    reg_alpha = 0.7, reg_lambda = 0.3,
                    subsample = 0.9, seed = 2401, print_every_n = 10,
                    silent = 1, nrounds = 800)#, watchlist=list(train=dtrain, test=dtest))

##%######################################################%##
#                                                          #
####                     read model                     ####
#                                                          #
##%######################################################%##


xgb.pred <- predict(xgb.best, dtrain)
#write_csv(as.data.frame(xgb.pred), "train_pred.csv")
postResample(xgb.pred, MeTrain$Active75)

confusionMatrix(table(as.numeric(xgb.pred > 0.5), MeTrain$Active75))

classifierplots(MeTrain$Active75, xgb.pred)

xgb.pred2 <- predict(xgb.best, dtest)
postResample(xgb.pred2, MeTest$Active75)

confusionMatrix(table(as.numeric(xgb.pred2>0.5), MeTest$Active75))

mat <- xgb.importance (feature_names = colnames(MeTrain),model = xgb.best)

xgb.plot.importance (importance_matrix = mat, top_n = 10, measure = NULL, rel_to_first = T,
                     main = "Relative Cell Morph Var Importance for the Model (Top 10)") 
write_csv(mat,"feature_importance_xgboost_active75_model.csv")

#selected best 30 from XGBOOST model

matop <- mat %>% arrange(desc(abs(Cover))) %>% top_n(10)
mabot <- mat %>% arrange(Cover) %>% top_n(-10)
bestCGP30 <- rbind(matop, mabot)
bestCGP30

############################################################
#                                                          #
#                       RPART Model                        #
#                                                          #
############################################################

library(rpart)
library(rpart.plot)

fit <- rpart(Active75 ~ .,
             method="class", data=ALL_valid_num_df[,-1])


printcp(fit) # display the results
plotcp(fit) # visualize cross-validation results
summary(fit) # detailed summary of splits

plot(fit, uniform=TRUE,
     main="Classification Tree for Active75")
text(fit, use.n=TRUE, all=TRUE, cex=.8)

# create attractive postscript plot of tree
post(fit, file = "modeltree.ps",
     title = "Classification Tree for Active75")


rpart.plot(fit)

############################################################
#                                                          #
#              plot some important parameters              #
#                                                          #
############################################################

library(reshape2)
library(ggforce)

ALL_valid_num_df_plot <- ALL_valid_num_df
ALL_valid_num_df_plot$Active75 <- if_else(ALL_valid_num_df_plot$Active75 == 1, "Active75", "Inactive75")

ALL_valid_num_df_plotm <- melt(ALL_valid_num_df_plot)

names(ALL_valid_num_df_plot)

ALL_valid_num_df_plotm %>% filter(variable != "abs_ID") %>% 
  filter(variable == "IntensityCellDigiPhaseContrastStdDev(DPC)") %>% 
  ggplot(., aes(Active75, value)) +
  geom_boxplot(aes(fill = Active75))

ALL_valid_num_df_plotm %>% filter(variable != "abs_ID") %>% 
  ggplot(., aes(Active75, value)) +
  geom_violin(aes(fill = Active75),
              draw_quantiles = 1,
              trim = TRUE,
              scale = "area") +
  #facet_wrap_(. ~ variable, ncol = 2, nrow = 2)
  facet_wrap_paginate(. ~ variable, nrow = 2, ncol = 2,
                      scales = "free", shrink = TRUE, labeller = "label_value",
                      as.table = TRUE, switch = NULL, drop = TRUE, dir = "h",
                      strip.position = "top", page = 1)

ALL_valid_num_df_plotm %>% filter(variable != "abs_ID") %>% 
  ggplot(., aes(Active75, value)) +
  geom_boxplot(aes(fill = Active75), outlier.shape = NA) +
  #facet_wrap_(. ~ variable, ncol = 2, nrow = 2)
  facet_wrap_paginate(. ~ variable, nrow = 2, ncol = 2,
                      scales = "free", shrink = TRUE, labeller = "label_value",
                      as.table = TRUE, switch = NULL, drop = TRUE, dir = "h",
                      strip.position = "top", page = 1)

ALL_valid_num_df_plotm %>% filter(variable != "abs_ID") %>% 
  ggplot(., aes(Active75, value)) +
  geom_boxplot(aes(fill = Active75), outlier.shape = NA) +
  #facet_wrap_(. ~ variable, ncol = 2, nrow = 2)
  facet_wrap_paginate(. ~ variable, nrow = 2, ncol = 2,
                      scales = "free", shrink = TRUE, labeller = "label_value",
                      as.table = TRUE, switch = NULL, drop = TRUE, dir = "h",
                      strip.position = "top", page = 2)


ALL_valid_num_df_plotm %>% filter(variable != "abs_ID") %>% 
  ggplot(., aes(Active75, value)) +
  geom_boxplot(aes(fill = Active75), outlier.shape = NA) +
  #facet_wrap_(. ~ variable, ncol = 2, nrow = 2)
  facet_wrap_paginate(. ~ variable, nrow = 2, ncol = 2,
                      scales = "free", shrink = TRUE, labeller = "label_value",
                      as.table = TRUE, switch = NULL, drop = TRUE, dir = "h",
                      strip.position = "top", page = 3)

ALL_valid_num_df_plotm %>% filter(variable != "abs_ID") %>% 
  ggplot(., aes(Active75, value)) +
  geom_boxplot(aes(fill = Active75), outlier.shape = NA) +
  #facet_wrap_(. ~ variable, ncol = 2, nrow = 2)
  facet_wrap_paginate(. ~ variable, nrow = 2, ncol = 2,
                      scales = "free", shrink = TRUE, labeller = "label_value",
                      as.table = TRUE, switch = NULL, drop = TRUE, dir = "h",
                      strip.position = "top", page = 4)

ALL_valid_num_df_plotm %>% filter(variable != "abs_ID") %>% 
  ggplot(., aes(Active75, value)) +
  geom_boxplot(aes(fill = Active75), outlier.shape = NA) +
  #facet_wrap_(. ~ variable, ncol = 2, nrow = 2)
  facet_wrap_paginate(. ~ variable, nrow = 2, ncol = 2,
                      scales = "free", shrink = TRUE, labeller = "label_value",
                      as.table = TRUE, switch = NULL, drop = TRUE, dir = "h",
                      strip.position = "top", page = 5)

ALL_valid_num_df_plotm %>% filter(variable != "abs_ID") %>% 
  ggplot(., aes(Active75, value)) +
  geom_boxplot(aes(fill = Active75), outlier.shape = NA) +
  #facet_wrap_(. ~ variable, ncol = 2, nrow = 2)
  facet_wrap_paginate(. ~ variable, nrow = 2, ncol = 2,
                      scales = "free", shrink = TRUE, labeller = "label_value",
                      as.table = TRUE, switch = NULL, drop = TRUE, dir = "h",
                      strip.position = "top", page = 6)

ALL_valid_num_df_plotm %>% filter(variable != "abs_ID") %>% 
  ggplot(., aes(Active75, value)) +
  geom_boxplot(aes(fill = Active75), outlier.shape = NA) +
  #facet_wrap_(. ~ variable, ncol = 2, nrow = 2)
  facet_wrap_paginate(. ~ variable, nrow = 2, ncol = 2,
                      scales = "free", shrink = TRUE, labeller = "label_value",
                      as.table = TRUE, switch = NULL, drop = TRUE, dir = "h",
                      strip.position = "top", page = 7)

Important <- as.vector(unique(ALL_valid_num_df_plotm$variable)[c(28,23,22,20,18,14,16,10,11,7,8,4)])

pdf("Selected_Cell_parameters_persingle_point.pdf", 8, 5)

for(i in seq(1:3)){ print(

  ALL_valid_num_df_plotm %>% filter(variable != "abs_ID") %>% 
    filter(variable %in% Important) %>% 
    ggplot(., aes(Active75, value)) +
    geom_boxplot(aes(fill = Active75), outlier.shape = NA) +
    #facet_wrap_(. ~ variable, ncol = 2, nrow = 2)
    facet_wrap_paginate(. ~ variable, nrow = 2, ncol = 2,
                        scales = "free", shrink = TRUE, labeller = "label_value",
                        as.table = TRUE, switch = NULL, drop = TRUE, dir = "h",
                        strip.position = "top", page = i)

)}
dev.off()

ALL_valid_num_df_plotm %>% filter(variable != "abs_ID") %>% 
  filter(variable %in% Important) %>% 
  ggplot(., aes(Active75, value)) +
  geom_boxplot(aes(fill = Active75), outlier.shape = NA) +
  #facet_wrap_(. ~ variable, ncol = 2, nrow = 2)
  facet_wrap_paginate(. ~ variable, nrow = 2, ncol = 2,
                      scales = "free", shrink = TRUE, labeller = "label_value",
                      as.table = TRUE, switch = NULL, drop = TRUE, dir = "h",
                      strip.position = "top", page = 2)

ALL_valid_num_df_plotm %>% filter(variable != "abs_ID") %>% 
  filter(variable %in% Important) %>% 
  ggplot(., aes(Active75, value)) +
  geom_boxplot(aes(fill = Active75), outlier.shape = NA) +
  #facet_wrap_(. ~ variable, ncol = 2, nrow = 2)
  facet_wrap_paginate(. ~ variable, nrow = 2, ncol = 2,
                      scales = "free", shrink = TRUE, labeller = "label_value",
                      as.table = TRUE, switch = NULL, drop = TRUE, dir = "h",
                      strip.position = "top", page = 3)


############################################################
#                                                          #
#                      only with IC50                      #
#                                                          #
############################################################

names(ALL_valid_v3)

ALL_valid_v3_IC50  <-  ALL_valid_v3 %>% filter(!is.na(IC50_nonumber)) %>% 
  select(c("abs_ID", "CPD_ID","CompoundNameB","Class", "IC50-roundGUF...10","IC50-roundGUF...11","IC50_nonumber",
           tidyselect::vars_select(names(ALL_valid_num), matches('?[\\)]'))))

ALL_valid_v3_IC50$Active75 <- if_else(ALL_valid_v3_IC50$`%Inhibition(DPC)` < 75, "inactive75", "active75")
ALL_valid_v3_IC50$IC50_class <- if_else(ALL_valid_v3_IC50$IC50_nonumber == ">25", "IC50>25", "IC50<25")


ALL_valid_v3_IC50m <- melt(ALL_valid_v3_IC50)

pdf("Selected_Cell_parameters_perIC50.pdf", 8, 5)

for(i in seq(1:12)){
  print(

  ALL_valid_v3_IC50m %>% filter(variable != "abs_ID") %>% filter(!is.na(Class)) %>% 
    filter(Class != "CYP") %>% 
    filter(Class != "PDE") %>% 
    filter(variable %in% Important) %>% 
    ggplot(., aes(IC50_class, value)) +
    geom_boxplot(aes(fill = IC50_class), outlier.shape = NA) +
    facet_wrap_paginate(variable ~ Class, nrow = 3, ncol = 3,
                        scales = "free", shrink = TRUE, labeller = "label_value",
                        as.table = TRUE, switch = NULL, drop = TRUE, dir = "h",
                        strip.position = "top", page = i)
  )}
dev.off()

############################################################
#                                                          #
#                   IC50 classification                    #
#                                                          #
############################################################

str(ALL_valid_v3_IC50)

ALL_valid_IC50_class <- ALL_valid_v3_IC50 %>% select(-c(abs_ID, Class, CompoundNameB, CPD_ID,
                                `IC50-roundGUF...10`, `IC50-roundGUF...11`, IC50_nonumber, 
                                Active75)) %>% select(-c(1,2)) %>% select(IC50_class, everything())


fitIC50 <- rpart(IC50_class ~ .,
             method="class", data=ALL_valid_IC50_class)


printcp(fitIC50) # display the results
plotcp(fitIC50) # visualize cross-validation results
summary(fitIC50) # detailed summary of splits

plot(fitIC50, uniform=TRUE,
     main="Classification Tree for IC_class")
text(fit, use.n=TRUE, all=TRUE, cex=.8)

# create attractive postscript plot of tree
post(fitIC50, file = "modeltree_IC50.ps",
     title = "Classification Tree for IC50_class")


rpart.plot(fitIC50)

############################################################
#                                                          #
#                         xgboost                          #
#                                                          #
############################################################

set.seed(2401)

ALL_valid_IC50_class$IC50_class <- if_else(ALL_valid_IC50_class$IC50_class == "IC50>25", 1,0)

trainIndex <- as.vector(createDataPartition(ALL_valid_IC50_class$IC50_class, p = .80, 
                                  list = FALSE, 
                                  times = 1))

MeTrain <- ALL_valid_IC50_class[ trainIndex,]
MeTest  <- ALL_valid_IC50_class[-trainIndex,]
sum(is.na(MeTrain[,-c(1)]))

labeltrain <- MeTrain$IC50_class
labeltest <- MeTest$IC50_class
table(labeltest)
table(labeltrain)

dtrain <- xgb.DMatrix(data = as.matrix(MeTrain[,-1]), label = labeltrain) 
dtest <- xgb.DMatrix(data = as.matrix(MeTest[,-1]),label = labeltest)

xgb.fit <- xgb.train(data = dtrain, label = labeltrain, booster = "gbtree",
                     objective = "binary:logistic",
                     colsample_bytree = 1, gamma = 0,
                     learning_rate = 0.01, max_depth = 6,
                     min_child_weight = 1,
                     reg_alpha = 0.8, reg_lambda = 0.2,
                     subsample = 0.8, print_every_n = 10,
                     silent = 1, nrounds = 500, watchlist=list(train=dtrain, test=dtest))

xgb.fit

cv.xgb <- xgb.cv(data = dtrain, nrounds = 190, nfold = 10, label = labeltrain,
                 objective = "binary:logistic", #booster = "gblinear",
                 colsample_bytree = 1, gamma = 0, booster = "gbtree",
                 learning_rate = 0.03, max_depth = 6,
                 min_child_weight = 1, #n_estimators = 7300,
                 reg_alpha = 0.8, reg_lambda = .2, print_every_n = 1,
                 subsample = 0.9)
cv.xgb

xgb.best <- xgboost(data = dtrain, label = labeltrain,
                    booster = "gbtree", objective = "binary:logistic",
                    colsample_bytree = 1, gamma = 0,
                    learning_rate = 0.01, max_depth = 6,
                    min_child_weight = 1, #n_estimators = 7300,
                    reg_alpha = 0.8, reg_lambda = 0.2,
                    subsample = 0.8, seed = 2401, print_every_n = 10,
                    silent = 1, nrounds = 500)#, watchlist=list(train=dtrain, test=dtest))

##%######################################################%##
#                                                          #
####                     read model                     ####
#                                                          #
##%######################################################%##


xgb.pred <- predict(xgb.best, dtrain)
postResample(xgb.pred, MeTrain$IC50_class)

confusionMatrix(table(as.numeric(xgb.pred > 0.5), MeTrain$IC50_class))

classifierplots(as.numeric(MeTrain$IC50_class), xgb.pred)
classifierplots_folder(MeTrain$IC50_class, xgb.pred, 
                       "W:/projects/ESP Projects Active/Corona_Broad_2020/GeneratedData/DataForDistribution/Figures_&_Tables")

xgb.pred2 <- predict(xgb.best, dtest)
postResample(xgb.pred2, MeTest$IC50_class)

confusionMatrix(table(as.numeric(xgb.pred2>0.5), MeTest$IC50_class))

mat <- xgb.importance (feature_names = colnames(MeTrain),model = xgb.best)

xgb.plot.importance (importance_matrix = mat, top_n = 10, measure = NULL, rel_to_first = T,
                     main = "Relative Cell Morph Var Importance for the Model (Top 10)") 
write_csv(mat,"feature_importance_xgboost_IC50_class_model.csv")

# PCA

library(ggrepel)
library(ggbiplot)

ALL_valid_IC50_class <- ALL_valid_IC50_class %>% select(-`NumberofAnalyzedFields(NucleiStaining)`)
ALL_valid_IC50_class$Compound <- ALL_valid_v3_IC50$CompoundNameB
ALL_valid_IC50_class$Class <- ALL_valid_v3_IC50$Class
ALL_valid_IC50_class$IC50_nonumber <- ALL_valid_v3_IC50$IC50_nonumber

PCA <- prcomp(ALL_valid_IC50_class[,-c(1,40:42)], scale = T,center = T) 
summary(PCA)

# check plot
screeplot(PCA, type = "l", npcs = 15, main = "Screeplot of the first 10 PCs")
abline(h = 1, col="red", lty=5)
legend("topright", legend=c("Eigenvalue = 1"),
       col=c("red"), lty=5, cex=0.6)

cumpro <- cumsum(PCA$sdev^2 / sum(PCA$sdev^2))
plot(cumpro[0:15], xlab = "PC #", ylab = "Amount of explained variance", main = "Cumulative variance plot")
abline(v = 6, col="blue", lty=5)
abline(h = 0.86911, col="blue", lty=5)
legend("topleft", legend=c("Cut-off @ PC6"),
       col=c("blue"), lty=5, cex=0.6)


library("factoextra")
library(plotly)
p <- fviz_pca_ind(PCA, geom.ind = "point", pointshape = 21, 
                  pointsize = 2, 
                  fill.ind = as.factor(ALL_valid_IC50_class$IC50_class), 
                  col.ind = "black", 
                  palette = "jco", 
                  addEllipses = F,
                  label = paste(ALL_valid_IC50_class$Compound, ALL_valid_IC50_class$IC50_nonumber, sep = " "),
                  col.var = "black",
                  repel = TRUE,
                  legend.title = "IC50") +
  ggtitle("2D PCA-plot from 38 cellular feature dataset") +
  theme(plot.title = element_text(hjust = 0.5))

ggplotly(p, tooltip = c("x","y","label") )

ALL_valid_v3_IC50$PC1 <- PCA$x[,1]
ALL_valid_v3_IC50$PC2 <- PCA$x[,2]

q <- ggplot(ALL_valid_v3_IC50, aes(PC1, PC2)) +
  geom_point(size = 1) +
  geom_point(aes(x = PC1, y = PC2, fill = as.factor(Class),
                 color = as.factor(Class))) +
  guides(size = F, fill = F) + labs(color = "Class") +
  theme(legend.position = c(0.9, 0.3),
        legend.direction = "vertical") +
  geom_text_repel(data=subset(ALL_valid_v3_IC50, IC50_class != "IC50<25"),
            aes(PC1, PC2, label=CompoundNameB))

ggplotly(q)

############################################################
#                                                          #
#                          purrr                           #
#                                                          #
############################################################


Final_total_annotated_original <- read_excel("202004_09_Final64_v2_16-19_35_59_total_annotated.xlsx", 
                                                         sheet = "Original")

Final_total_annotated_original_2merge <- Final_total_annotated_original %>% select(`%Inhibition Cells...Number.of.Objects_norm`,
                                                                                       `Compound ID`, `Customer ID`,
                                                                                      `Code_Activity`, Name, `NAME from BROAD`)

sum(is.na(Final_total_annotated_original_2merge$`Customer ID`))

names(Final_total_annotated_original_2merge)

names(Final_total_annotated_original_2merge)[2] <- "CompoundId"
names(Final_total_annotated_original_2merge)[3] <- "CPDID"

ALL_valid_v3_purrr <- inner_join(ALL_valid_v3, Final_total_annotated_original_2merge, by = "CompoundId")
ALL_valid_v3_purrr_num <- ALL_valid_v3_purrr %>% select_if(is.numeric) %>% select(-c(1:4, 46))

ALL_valid_v3_purrr_num$Code_Activity <- ALL_valid_v3_purrr$Code_Activity

ALL_valid_v3_purrr_num <- ALL_valid_v3_purrr_num %>% select(Code_Activity, everything())

unique(ALL_valid_v3_purrr_num$Code_Activity)

ALL_valid_v3_purrr_num$Code <- if_else(ALL_valid_v3_purrr_num$Code_Activity == "Active_75", 1, 
                                       if_else(ALL_valid_v3_purrr_num$Code_Activity == "Active_50", 2,
                                               if_else(ALL_valid_v3_purrr_num$Code_Activity == "Inactive_25", 3,
                                                       if_else(ALL_valid_v3_purrr_num$Code_Activity == "Inactive_0", 4,
                                                               if_else(ALL_valid_v3_purrr_num$Code_Activity == "Inactive_neg", 5,6)))))


unique(ALL_valid_v3_purrr_num$Code)
table(ALL_valid_v3_purrr_num$Code)

ALL_valid_v3_purrr_num$Code_Activity <- NULL

ALL_valid_v3_purrr_num %>% 
  map(~lm(ALL_valid_v3_purrr_num$Code ~ .x, data = ALL_valid_v3_purrr_num)) %>% 
  map(summary) %>% 
  map_dbl("r.squared") %>% 
  tidy %>% 
  dplyr::arrange(desc(x)) %>% 
  rename(r.squared = x) -> r2s

kable(r2s)

#### Plot melted 

library(reshape2)
library(ggforce)

ALL_valid_v3_purrr_num$Code <- as.factor(ALL_valid_v3_purrr_num$Code)

ALL_valid_v3_purrr_num_m <- melt(ALL_valid_v3_purrr_num, )

names(ALL_valid_v3_purrr_num_m)

Important <- names(ALL_valid_v3_purrr_num)[1:41]

ALL_valid_v3_purrr_num_m %>% filter(variable != "abs_ID") %>% 
  filter(variable == "NumberOfNormalNuclei(NucleiStaining)") %>% 
  ggplot(., aes(Code, value)) +
  geom_boxplot(aes(fill = Code)) +
  scale_x_discrete(labels=c("1" = ">75% inh", "2" = "50< %inh <=75","3" = "25< %inh <50", "4" = "0< %inh <25", 
                          "5" = "-10< %inh < 0", "6" = "%inh < -10"))

pdf("Selected_Cell_parameters_per_CODE.pdf", 8, 5)

for(i in seq(1:41)){
  print(
    ALL_valid_v3_purrr_num_m %>% 
      filter(variable %in% Important) %>% filter(variable != "%Confluence(DPC)") %>% 
      ggplot(., aes(Code, value)) +
      geom_boxplot(aes(fill = Code), outlier.shape = NA) +
      scale_x_discrete(labels=c("1" = ">75% inh", "2" = "50< %inh <=75","3" = "25< %inh <50", "4" = "0< %inh <25", 
                                "5" = "-10< %inh < 0", "6" = "%inh < -10")) +
      facet_wrap_paginate(variable ~ ., nrow = 2, ncol = 1,
                          scales = "free", shrink = TRUE, labeller = "label_value",
                          as.table = TRUE, switch = NULL, drop = TRUE, dir = "h",
                          strip.position = "top", page = i)
  )}


dev.off()

pdf("Selected_Cell_parameters_per_CODE_nopotentiator.pdf", 8, 5)

for(i in seq(1:41)){
  print(
    ALL_valid_v3_purrr_num_m %>% 
      filter(variable %in% Important) %>% filter(variable != "%Confluence(DPC)") %>% 
      filter(Code != "6") %>% droplevels() %>% 
      ggplot(., aes(Code, value)) +
      geom_boxplot(aes(fill = Code), outlier.shape = NA) +
      scale_x_discrete(labels=c("1" = ">75% inh", "2" = "50< %inh <=75","3" = "25< %inh <50", "4" = "0< %inh <25", 
                                "5" = "-10< %inh < 0", "6" = "%inh < -10")) +
      facet_wrap_paginate(variable ~ ., nrow = 2, ncol = 1,
                          scales = "free", shrink = TRUE, labeller = "label_value",
                          as.table = TRUE, switch = NULL, drop = TRUE, dir = "h",
                          strip.position = "top", page = i)
  )}


dev.off()


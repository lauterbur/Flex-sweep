library(tidyverse)
library(yardstick)
library(viridis)

args = commandArgs(trailingOnly=TRUE)
neutral_preds <- args[1]
sweep_preds <- args[2]
outputDir <- args[3]

### load data
neutral_data <- read_tsv(neutral_preds)
sweep_data <- read_tsv(sweep_preds)
pred_data <- bind_rows(neutral = neutral_data, sweep = sweep_data, .id = "type")

print(neutral_preds)
### get true/false counts
confusion_data <- pred_data %>%
  group_by(type) %>%
  count(predicted_class,.drop=FALSE) %>%
  mutate(true_false=ifelse(type==predicted_class & type=="neutral","true_negative",
                           ifelse(type==predicted_class & type=="sweep","true_positive",
                                  ifelse(type!=predicted_class & type=="neutral","false_positive","false_negative")))) %>%
  ungroup(type) %>%
  select(-c(type,predicted_class)) %>%
  spread(true_false,n,fill=0,drop=FALSE) #%>%

## get rates
rate_data <- confusion_data %>%
  ungroup() %>%
  mutate(false_negative=false_negative/rowSums(.[c(1,4)]),
         false_positive=false_positive/rowSums(.[c(2,3)]),
         true_negative=true_negative/rowSums(.[c(2,3)]),
         true_positive=true_positive/rowSums(.[c(1,4)])) %>%
  mutate(false_negative=replace_na(false_negative, 0),
         false_positive=replace_na(false_positive, 0),
         true_positive=replace_na(true_positive, 0),
         true_negative=replace_na(true_negative, 0)) %>%
  mutate(accuracy=(true_positive+true_negative)/(true_positive+true_negative+false_positive+false_negative),
         precision=true_positive/(true_positive + false_positive))

## set up data for ROC
roc_data <- pred_data %>% 
  mutate(type=relevel(factor(type),"sweep")) %>% # because the first level needs to be the event of interest
  roc_auc(type,`prob(sweep)`) %>%
  rename(AUC=`.estimate`) %>%
  select(-c(`.metric`,`.estimator`))

## join data
rate_roc_data <- bind_cols(rate_data,roc_data)
pred_rate_auc_data<-bind_cols(pred_data,rate_roc_data)

write_csv(rate_roc_data, paste0(outputDir,"/plots/training_classification_rates.csv"))

p <- pred_rate_auc_data %>%
  mutate(type=relevel(factor(type),"sweep")) %>% # because the first level needs to be the event of interest
  group_modify(~roc_curve(.x,type,`prob(sweep)`),.keep=TRUE) %>%
  ggplot(.) + 
  aes(x=1-specificity, y=sensitivity) +
  geom_path(size=2, color="orange") +
  coord_equal() + 
  theme_bw() +
  labs(x="false positive rate",y="power")

fprs <- pred_rate_auc_data %>%
  mutate(type=relevel(factor(type),"sweep")) %>%
  select(false_positive,false_negative,accuracy,precision) %>% 
  distinct() %>%
  mutate_if(is.numeric,round,4)

p <- p +
  geom_text(x=.2,y=.6,aes(label="FPR"), color="black") +
  geom_text(x=.2,y=.5,aes(label=round(false_positive,4)),data=fprs) +
  geom_text(x=.43,y=.6,aes(label="FNR"), color="black") +
  geom_text(x=.43,y=.5,aes(label=round(false_negative,4)),data=fprs) +
  geom_text(x=.66,y=.6,aes(label="Accuracy"), color="black") +
  geom_text(x=.66,y=.5,aes(label=round(accuracy,4)),data=fprs) +
  geom_text(x=.9,y=.6,aes(label="Precision"), color="black") +
  geom_text(x=.9,y=.5,aes(label=round(precision,4)),data=fprs) 
  
ggsave(paste0(outputDir,"/plots/training_ROC.png"), p)

### load history data
history <- paste0(outputDir,"/",outputDir,"Model/",outputDir,"_history.csv")
history_data <- read_csv(history)

h <- history_data %>%
  select(loss, val_loss, accuracy, val_accuracy) %>%
  mutate(epoch=row_number()) %>%
  pivot_longer(cols=c(loss, val_loss, accuracy, val_accuracy), names_to="metric_name", values_to="metric_val") %>%
  separate(metric_name, c("A","B")) %>%
  mutate(B=ifelse(is.na(B), A, B),
         A=ifelse(A!="val", "train", A)) %>%
  group_by(A, B) %>%
  ggplot() + aes(x=epoch, y = metric_val) +
  geom_line(aes(color=B, lty=A), size=1) +
  theme_bw() +
  scale_color_manual(name = "", values=c("blue","orange")) +
  scale_linetype_manual(name = "", values=c("solid", "dotted"), labels=c("train", "validation")) +
  ylab("value") +
  labs(title = "history")

ggsave(paste0(outputDir,"/plots/training_history.png"), h)

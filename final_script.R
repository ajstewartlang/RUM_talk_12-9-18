# First load the tidyverse packages
library(tidyverse)

# Load the packages we need for our linear mixed models
library(lme4)
library(lmerTest)
library(emmeans)

# Load the data 
data <- read_csv("fpfp0.ixs")

# Labelling for each condition is as follows:
# C1 = Neutral condition
# C2 = Negative condition
# C3 = Positive condition

# Relabel the Condition variable ("cond") and set as factor
data$cond <- recode (data$cond, "C1"="Neutral", "C2"="Negative", "C3"="Positive")
data$cond <- as.factor(data$cond)

# Filter so we only have the Region 3 data - this is our critical analysis region
# Ignore any zero reading times - these will be tracker loss/missing data
data_R3 <- filter (data, reg=="R3" & DV > 0)

data_R3

# Let's visualise the data for each of the 3 conditions
ggplot (data_R3, aes(x=DV)) + geom_histogram() + facet_wrap(~cond)

# Let's look at the data on a Participant by Participant basis
ggplot (data_R3, aes (x=cond, y=DV, colour=cond)) + geom_jitter(width=.1, alpha=.5) + 
  stat_summary(fun.y=mean, geom="point", colour="black") +
  facet_wrap(~subj) +
  scale_color_discrete(guide=FALSE) + 
  labs (y="Reading Time (msec.)", x="Condition", 
        title="Critical Region, First Pass Reading Times Facted Wrapped by Participant")

# Let's work out for how many participants the effect went in the direction as predicted
# First we're aggregating by participants
data_agg <- data_R3 %>% group_by(subj, cond) %>% summarise (mean=mean(DV), sd=sd(DV))

# For how many people was the Negative condition faster than the Positive?
sum (data_agg[data_agg$cond=="Negative",]$mean < data_agg[data_agg$cond=="Positive",]$mean)

# For how many people was the Negative condition faster than the Neutral?
sum (data_agg[data_agg$cond=="Negative",]$mean < data_agg[data_agg$cond=="Neutral",]$mean)

# For how many people were the above both TRUE? 
sum (data_agg[data_agg$cond=="Negative",]$mean < data_agg[data_agg$cond=="Positive",]$mean & 
       data_agg[data_agg$cond=="Negative",]$mean < data_agg[data_agg$cond=="Neutral",]$mean)

# Let's do the same but this time by Items 
ggplot (data_R3, aes (x=cond, y=DV, colour=cond)) + geom_jitter(width=.1, alpha=.5) + 
  stat_summary(fun.y=mean, geom="point", colour="black") +
  facet_wrap(~item) +
  scale_color_discrete(guide=FALSE) + 
  labs (y="Reading Time (msec.)", x="Condition",
                title="Critical Region, First Pass Reading Times Facted Wrapped by Item")

data_agg <- data_R3 %>% group_by(item, cond) %>% summarise (mean=mean(DV), sd=sd(DV))

sum (data_agg[data_agg$cond=="Negative",]$mean < data_agg[data_agg$cond=="Positive",]$mean)
sum (data_agg[data_agg$cond=="Negative",]$mean < data_agg[data_agg$cond=="Neutral",]$mean)

# Let's build a raincloud plot
library(RColorBrewer)
library(plyr) #note, need to detach this after this plot as clashes with aspects of dplyr
source("https://gist.githubusercontent.com/ajstewartlang/6c4cd8ab9e0c27747424acdfb3b4cff6/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

raincloud_theme = theme(
  text = element_text(size = 12),
  axis.title.x = element_text(size = 12),
  axis.title.y = element_text(size = 12),
  axis.text = element_text(size = 12),
  axis.text.x = element_text(angle = 45, vjust = 0.5),
  legend.title=element_text(size=12),
  legend.text=element_text(size=12),
  legend.position = "right",
  #plot.title = element_text(lineheight=.8, size = 12),
  panel.border = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))

lb <- function(x) mean(x) - sd(x)
ub <- function(x) mean(x) + sd(x)

sumld <- ddply(data_R3, ~cond, summarise, mean = mean(DV), median = median(DV), lower = lb(DV), upper = ub(DV))
head(sumld)

ggplot(data = data_R3, aes(y = DV, x = cond, fill = cond)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8, trim=FALSE) +
  geom_point(aes(y = DV, color = cond), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = .1,  outlier.shape = NA, alpha = 0.5) +
  #expand_limits(x = 3) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  scale_color_brewer(palette = "Accent") +
  scale_fill_brewer(palette = "Accent") +
  coord_flip() +
  theme_bw() +
  raincloud_theme +
  labs(y="Reading Time (msec.)", x="Condition", title="Critical Region, First Pass Reading Times")

detach("package:plyr", unload=TRUE) #need to detach as some clashes with dplyr

# Build our first linear mixed model with maximal random effects structure
model.null <- lmer(DV ~ (1+cond|subj) + (1+cond|item), data_R3, REML=TRUE)
model <- lmer(DV ~ cond + (1+cond|subj) + (1+cond|item), data_R3, REML=TRUE)

# Check the experimental model differs from the Null model
anova (model.null, model)

# Generate a summary of our model parameters plus run some pairwise comparisons 
summary (model)
emmeans (model, pairwise~cond, adjust="none")

# Check we're happy witht the model residuals
qqnorm (residuals(model))

# Now let's look at the subsequent region of text (Region 4)
# Exclude missing data/tracker loss points where the DV==0
data_R4 <- filter (data, reg=="R4" & DV > 0)

# Need to load the plyr package again
library(plyr)

sumld <- ddply(data_R4, ~cond, summarise, mean = mean(DV), median = median(DV), lower = lb(DV), upper = ub(DV))
head(sumld)

ggplot(data = data_R4, aes(y = DV, x = cond, fill = cond)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8, trim=FALSE) +
  geom_point(aes(y = DV, color = cond), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = .1,  outlier.shape = NA, alpha = 0.5) +
  #expand_limits(x = 3) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  scale_color_brewer(palette = "Accent") +
  scale_fill_brewer(palette = "Accent") +
  coord_flip() +
  theme_bw() +
  raincloud_theme +
  labs(y="Reading Time (msec.)", x="Condition", title="Post-Critical Region, First Pass Reading Times")

detach("package:plyr", unload=TRUE) #need to detach as some clashes with dplyr

# Build the linear mixed model for Region 4 as we did for Region 3
model.null <- lmer(DV ~ (1|subj) + (1|item), data_R4, REML=TRUE)
model <- lmer(DV ~ cond + (1|subj) + (1|item), data_R4, REML=TRUE)

anova (model.null, model)

summary (model)
emmeans (model, pairwise~cond, adjust="none")

# Read in the Regressions Out data
FPRO <- read_csv("fproro0.ixs")

# Recode as before
data_R3 <- filter (FPRO, reg=="R3" & cond != "C0")
data_R3$cond <- recode (data_R3$cond, "C1"="Neutral", "C2"="Negative", "C3"="Positive")
data_R3$cond <- as.factor (data_R3$cond)

# Let's generate some summary data
data_summ <- data_R3 %>% group_by(cond) %>% summarise(mean=mean(DV), sd=sd(DV))

# Build a bar chart
ggplot (data_summ, aes (x=cond, y=mean, fill=cond)) + geom_col() + 
  geom_text(aes (label= format(round(mean, 0))), nudge_y = 1) +
  scale_fill_discrete (guide=FALSE) + labs (y="Regressions (%)", x="Condition", 
                                            title="Regressions Out of Critical Region")

# Recode the DV as binomial
data_R3$DV <- recode (data_R3$DV, "100"=1 ,"0"=0)

# Build the null and experimental glmm given sampling from the binomial distribution
model.null <- glmer (DV ~ (1|subj) + (1|item), data_R3, family=binomial)
model <- glmer (DV ~ cond + (1|subj) + (1|item), data_R3, family=binomial)
anova (model.null, model)
summary (model)

# Now look at Regressions Out of the post-critical region (Region 4)
data_R4 <- filter (FPRO, reg=="R4" & cond != "C0")
data_R4$cond <- recode (data_R4$cond, "C1"="Neutral", "C2"="Negative", "C3"="Positive")
data_R4$cond <- as.factor (data_R4$cond)

data_summ <- data_R4 %>% group_by(cond) %>% summarise(mean=mean(DV), sd=sd(DV))

# Build a bar chart
ggplot (data_summ, aes (x=cond, y=mean, fill=cond)) + geom_col() + 
  geom_text(aes (label= format(round(mean, 0))), nudge_y = 1) +
  scale_fill_discrete (guide=FALSE) + labs (y="Regressions (%)", x="Condition",
                                            title="Regressions Out of Post-Critical Region")

# Recode DV as binomial
data_R4$DV <- recode (data_R4$DV, "100"=1 ,"0"=0)

# Build the glmms
model.null <- glmer (DV ~ (1|subj) + (1|item), data_R4, family=binomial)
model <- glmer (DV ~ cond + (1|subj) + (1|item), data_R4, family=binomial)

anova (model.null, model)
summary (model)

# As model is significant, run some pairwise comparisions to determine where the effect is. Ask for the
# output on the response scale.
emmeans (model, pairwise~cond, adjust="none", type="response")

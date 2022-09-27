# R code used in:
# Towards reproducible analysis of benthos structural complexity: 
# A case study on Antarctic polychaete reefs using action cameras

# Authors: J.C. Montes-Herrera, G. Johnstone, J. Stark, N. Hill, V. Cummings, V. Lucieer

# Contact: juancarlos.montesherrera@utas.edu.au

###### Linear mixed models ######

library(ggplot2)
library(lme4)

path<- "C:\\Users\\jcmontes\\Desktop\\"

# Load file:
df<- read.csv(paste0(path, "habitat_metrics.csv"))

# ----

# Fractal dimension (D32) linear models:
d32_pc_lm<- lm(fd32 ~ Polychaete_Colonies, data=df)
d32_mud_lm<- lm(fd32 ~ Mud, data=df)
d32_bpt_lm<- lm(fd32 ~ Broken_Polychaete_Tubes, data=df)

# Profile Curvature linear models:
profc_pc_lm<- lm(mean_profile_curvature ~ Polychaete_Colonies, data=df)
profc_mud_lm<- lm(mean_profile_curvature ~ Mud, data=df)
profc_bpt_lm<- lm(mean_profile_curvature ~ Broken_Polychaete_Tubes, data=df)

# ----

# D32 linear mixed model with Transect as a random effect:
d32_pc_lme <- lmer(fd32 ~ Polychaete_Colonies + (1 | Transect),
                   data = df)
d32_mud_lme <- lmer(fd32 ~ Mud + (1 | Transect),
                    data = df)
d32_bpt_lme <- lmer(fd32 ~ Broken_Polychaete_Tubes + (1 | Transect),
                    data = df)

# Profile curvature linear mixed model with Transect as a random effect:
profc_pc_lme <- lmer(mean_profile_curvature ~ Polychaete_Colonies + (1 | Transect),
                     data = df)
profc_mud_lme <- lmer(mean_profile_curvature ~ Mud + (1 | Transect),
                      data = df)
profc_bpt_lme <- lmer(mean_profile_curvature ~ Broken_Polychaete_Tubes + (1 | Transect),
                      data = df)

# ----

# D32 linear mixed model with Random intercept and slope

d32_pc_rslme <- lmer(fd32 ~ Polychaete_Colonies + (1 + Polychaete_Colonies|Transect), 
                      data = df)

d32_mud_rslme <- lmer(fd32 ~ Mud + (1 + Mud|Transect), 
                      data = df)

d32_bpt_rslme <- lmer(fd32 ~ Broken_Polychaete_Tubes + (1 + Broken_Polychaete_Tubes|Transect), 
                       data = df)

# Profile Curvature linear mixed model with Random intercept and slope

profc_pc_rslme <- lmer(mean_profile_curvature ~ Polychaete_Colonies + 
                        (1 + Polychaete_Colonies |Transect), 
                      data = df)

profc_mud_rslme <- lmer(mean_profile_curvature ~ Mud + 
                        (1 + Mud |Transect), 
                      data = df)

profc_bpt_rslme <- lmer(mean_profile_curvature ~ Broken_Polychaete_Tubes + 
                        (1 + Broken_Polychaete_Tubes|Transect), 
                      data = df)

# ----
# Model decision based on AIC

anova(profc_bpt_rslme, profc_bpt_lme, profc_bpt_lm)

## Performance and goodness-of-fit
performance::r2(d32_pc_lme)

# Plots

library(ggeffects)
pred_mm <- ggpredict(profc_bpt_lme, terms = c("Broken_Polychaete_Tubes")) # extract prediction dataframe

p <- (ggplot(pred_mm) + 
    geom_line(aes(x = x, y = predicted), color="black", size=1) +          # slope
    geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
                fill = "darkgrey", alpha = 0.5) +  # error band
    geom_point(data = df,                      # adding the raw data (scaled values)
               aes(x = Broken_Polychaete_Tubes, y = mean_profile_curvature, color = Transect), size=2) + 
    labs(x = "Broken_Polychaete_Tubes (% cover)", y = "Profile Curvature", 
         title = "Broken_Polychaete_Tubes - Model 2") + 
    theme(legend.box.background = element_blank())
)


p

#Save figure
tiff(file="C:\\Users\\jcmontes\\Desktop\\profc-bpt_M2.png", 
     units="in", width=5, height=4, res=300)
p
dev.off()



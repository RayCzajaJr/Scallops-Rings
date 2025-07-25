
library(dplyr)
library(lubridate)
library(MuMIn)
library(betareg)
library(ggridges)
library(lme4)
library(stats)
library(car)
library(viridis) 
library(DHARMa)
library(ggeffects)
library(effects)
library(mgcv)
library(DescTools)
library(rcartocolor)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)


####### OBJECTIVE: Compare mean/median sizes at sites where we have multiple yrs of data 


# Rename the DataFrame
sizedistscallops <- X_Master_ring_total_shell_ht_22Feb2025

unique_locations <- unique(sizedistscallops$Location)
print(unique_locations)

# Create the new 'location' column with grouped labels
sizedistscallops <- sizedistscallops %>%
  mutate(location = case_when(
    Location %in% c("East Marion", "East Marion (NE of Bay Ave)") ~ "East Marion",
    Location %in% c("Southold (Off Cedar Beach)", "Southold Bay (off Cedar Beach)", "Southold Bay (N central)") ~ "Southold - Cedar Beach",
    Location %in% c(
                    "NW Harbor (East side)", "NW Harbor (Off Split Rock)", 
                    "NW Harbor (S of Alewife Creek)", "NW Harbor (South of Alewife Creek)") ~ "NW Harbor - E Side",
    Location %in% c("NW Harbor (Barcelona Point)",  "NW Harbor (Barcelona Pt)", "NWHarbor (N of Mile Hill Road)") ~ "NW Harbor - Barcelona ",
    Location == "Flanders" ~ "Flanders",
    Location == "Noyack Bay (EW Side)" ~ "Noyack Bay",
    Location %in% c("Hallock Bay", "Hallock Bay (Outside Narrow River)", "Hallock Bay (central flats)") ~ "Hallock Bay",
    Location %in% c("Hog Neck", "Hog Neck Bay (SE corner)") ~ "Hog Neck",
    Location %in% c("OH Harbor", "OH Harbor (North)", "OH Harbor North", "OH Harbor - N") ~ "OH Harbor",
    Location %in% c("Robins Island (West central side)", "Robins Island (W side)") ~ "Robins Island",
    Location %in% c("Shelter Island (Hay Beach)") ~ "Shelter Island - H Beach",
    Location %in% c("Shelter Island (East Side)", "Shelter Island (NE Side)") ~ "Shelter Island - E Side",
    TRUE ~ Location # Keep original value if no match
  ))

unique_locations <- unique(sizedistscallops$location)
print(unique_locations)

sizedistscallops$sample <- paste(sizedistscallops$Location, sizedistscallops$date, sep = "_")

unique_sample <- unique(sizedistscallops$sample)
print(unique_sample)

# Convert date column to Date format
sizedistscallops <- sizedistscallops %>%
  mutate(date = mdy(date),  
         year = year(date)) 

# Make and apply Function to get the first full week of oct for a given year
get_first_full_week_oct <- function(year) {
  october_days <- seq(ymd(paste0(year, "-10-01")), ymd(paste0(year, "-10-31")), by = "day")
    first_monday <- october_days[wday(october_days) == 2][1]  
    if (is.na(first_monday)) return(as.Date(character()))
    week_range <- seq(first_monday, first_monday + 6, by = "day")
  return(week_range)
}

first_full_weeks_oct <- do.call(rbind, lapply(unique(sizedistscallops$year), function(y) {
  data.frame(year = y, date = get_first_full_week_oct(y))
}))

first_full_weeks_oct$date <- as.Date(first_full_weeks_oct$date) 

sizedistscallops_oct <- sizedistscallops %>%
  inner_join(first_full_weeks_oct, by = c("year", "date"))

# Make and apply Function to get the first full week of nov for a given year
get_first_full_week_nov <- function(year) {
  november_days <- seq(ymd(paste0(year, "-11-01")), ymd(paste0(year, "-11-30")), by = "day")
    first_monday <- november_days[wday(november_days) == 2][1]  
    if (is.na(first_monday)) return(as.Date(character()))
    week_range <- seq(first_monday, first_monday + 6, by = "day")
  return(week_range)
}

first_full_weeks_nov <- do.call(rbind, lapply(unique(sizedistscallops$year), function(y) {
  data.frame(year = y, date = get_first_full_week_nov(y))
}))

first_full_weeks_nov$date <- as.Date(first_full_weeks_nov$date) 

sizedistscallops_nov <- sizedistscallops %>%
  inner_join(first_full_weeks_nov, by = c("year", "date"))

# Convert year to character for analyses 
sizedistscallops_nov <- sizedistscallops_nov %>%
  mutate(Year = as.character(Year))

sizedistscallops_oct <- sizedistscallops_oct %>%
  mutate(Year = as.character(Year))

sizedistscallops<- sizedistscallops %>%
  mutate(Year = as.character(Year))

# Ridge plot for shell and ring heights by site and year (two seperte months for shell heights)
novtotalhtyear<-ggplot(sizedistscallops_nov, aes(x = `Total Ht`, y = Year, fill = Year)) +
  stat_density_ridges(quantile_lines = TRUE, alpha = 0.7,
                      quantiles = 2)+
  scale_fill_carto_d(palette = "BurgYl") + 
  theme_bw() +  
  labs(
       x = "Total Height (Nov)",
       y = "Year") +
  theme(legend.position = "none",  
        plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(limits = c(53, NA))

novtotalhtloc<-ggplot(sizedistscallops_nov, aes(x = `Total Ht`, y = location, fill = location)) +
  stat_density_ridges(quantile_lines = TRUE, alpha = 0.7,
                      quantiles = 2)+
  scale_fill_carto_d(palette = "BluGrn") + 
  theme_bw() +  
  labs(
       x = "Total Height (Nov)",
       y = "Site") +
  theme(legend.position = "none",  
        plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(limits = c(53, NA))

octtotalhtyear<-ggplot(sizedistscallops_oct, aes(x = `Total Ht`, y = Year, fill = Year)) +
  stat_density_ridges(quantile_lines = TRUE, alpha = 0.7,
                      quantiles = 2)+
  scale_fill_carto_d(palette = "BurgYl") + 
  theme_bw() +  
  labs(
       x = "Total Height (Oct)",
       y = "Year") +
  theme(legend.position = "none",  
        plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(limits = c(46, 88))

octtotalhtloc<-ggplot(sizedistscallops_oct, aes(x = `Total Ht`, y = location, fill = location)) +
  stat_density_ridges(quantile_lines = TRUE, alpha = 0.7,
                      quantiles = 2)+
  scale_fill_carto_d(palette = "BluGrn") + 
  theme_bw() +  
  labs(
       x = "Total Height (Oct)",
       y = "Site") +
  theme(legend.position = "none",  
        plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(limits = c(46, 88))

ridgeplots<-grid.arrange(octtotalhtloc, octtotalhtyear, novtotalhtloc, novtotalhtyear,
                                 ncol = 2, nrow = 2)

ggsave("ridgeplots.tiff",ridgeplots, dpi = 300, bg = "white",
       width = 24,
       height = 24,
       units = "cm")

ringhtloc<-ggplot(sizedistscallops, aes(x = `Ring Ht`, y = location, fill = location)) +
  stat_density_ridges(quantile_lines = TRUE, alpha = 0.7,
                      quantiles = 2)+
  scale_fill_carto_d(palette = "BluGrn") + 
  theme_bw() +  
  labs(
       x = "Ring Height",
       y = "Site") +
  theme(legend.position = "none",  
        plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(limits = c(0, 80))

# Boxplot since there are too many years for a nice ridge plot
ringhtyr<-ggplot(sizedistscallops, aes(x = as.factor(Year), y = `Ring Ht`, fill = as.factor(Year))) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +  
  scale_fill_carto_d(palette = "BurgYl") +  
  theme_bw() +  
  labs(
       x = "Year",  
       y = "Ring Height") +
  theme(legend.position = "none",  
        plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0, 78))

ringhts<-grid.arrange(ringhtloc, ringhtyr,
                         ncol = 1, nrow = 2)

ggsave("ringhts.tiff",ringhts, dpi = 300, bg = "white",
       width = 20,
       height = 24,
       units = "cm")

# Prelim analyses for differences in shell and ring ht by year and location.... may wanna switch to a KS test?
moct <- aov(sizedistscallops_oct$`Total Ht`~ location + Year + location * Year, data=sizedistscallops_oct)
summary(moct)
res1 <- simulateResiduals(moct)
plot(res1)
mnov<- aov(sizedistscallops_nov$`Total Ht`~ location + Year + location * Year, data=sizedistscallops_nov)
summary(mnov)
mring <- aov(sizedistscallops$`Ring Ht`~ location + Year + location * Year, data=sizedistscallops)
summary(mring)


####### OBJECTIVE: Examine relationship between growth after the winter vs ring size (ie size reached before the winter)


# Add a month column to each dataframe before merging
sizedistscallops_oct <- sizedistscallops_oct %>%
  mutate(Month = "October")

sizedistscallops_nov <- sizedistscallops_nov %>%
  mutate(Month = "November")

# Merge the two monthly dataframes into one
sizedistscallops_octandnov <- bind_rows(sizedistscallops_oct, sizedistscallops_nov)

colnames(sizedistscallops_octandnov) <- gsub(" ", "_", colnames(sizedistscallops_octandnov))

m <- lmer(Growthafterfirstwinter~ Ring_Ht + (1 | location) + (1 | Year), 
          data = sizedistscallops_octandnov)
Anova(m)
r.squaredGLMM(m)
res1 <- simulateResiduals(m)
plot(res1)

efct <- effect("Ring_Ht", mod = m, xlevels = list(Ring_Ht = seq(0, 72, length.out = 100)))

efct <- as.data.frame(efct)

growthpostwinter<-ggplot() +
  geom_point(data=sizedistscallops_octandnov, aes(x = Ring_Ht, y =Growthafterfirstwinter), alpha = 0.4, color = "black")+  
  geom_line(data = efct, aes(x = Ring_Ht, y = fit), color = "darkseagreen4", size=1) +
  geom_ribbon(data = efct, aes(x = Ring_Ht, ymin = lower, ymax = upper),  fill = "darkseagreen4", alpha = 0.3) +
  theme_bw() +
  labs(
       x = "Ring Height",
       y = "Growth after 1st Winter") +
  theme(plot.title = element_text(hjust = 0.5))  

ggsave("growthpostwinter.tiff",growthpostwinter, dpi = 300, bg = "white",
       width = 20,
       height = 15,
       units = "cm")


####### OBJECTIVE: Examine relationship between % ripe in Oct/Nov vs small ring adults the next yea


scallopringsummarydf <- X_2024_master_scallop_ring_paper_summary_data_highlighting_revised15Mar2025_data_UNCHANGED_4

# Create column for days post 9-30
scallopringsummarydf$Date <- as.Date(scallopringsummarydf$Date)
scallopringsummarydf$dayspostsep30 <- as.numeric(scallopringsummarydf$Date - as.Date(paste0(format(scallopringsummarydf$Date, "%Y"), "-09-30")))

# Convert percent to decimal and remove NAS for percent ripe
scallopringsummarydf <- scallopringsummarydf[!is.na(scallopringsummarydf$orangegonad_percentripe), ]
scallopringsummarydf$orangegonad_percentripe <- scallopringsummarydf$orangegonad_percentripe / 100

# Convert percent to decimal and remove NAS for percent small rings
scallopringsummarydf <- scallopringsummarydf[!is.na(scallopringsummarydf$percentsmallrings_lessthan20mm), ]
scallopringsummarydf$percentsmallrings_lessthan20mm <- scallopringsummarydf$percentsmallrings_lessthan20mm / 100

# Remove months after the fall
scallopringsummarydf <- scallopringsummarydf[scallopringsummarydf$dayspostsep30 <= 60, ]

# Get annnaul means with lagged predictor
scallopsummarydf_annualmeans <- scallopringsummarydf %>%
  mutate(year = format(Date, "%Y")) %>%
  group_by(year) %>%
  summarise(
    mean_orangegonad = mean(orangegonad_percentripe, na.rm = TRUE),
    mean_percentsmallrings = mean(percentsmallrings_lessthan20mm, na.rm = TRUE)
  ) %>%
  ungroup()

scallopsummarydf_annualmeans <- scallopsummarydf_annualmeans %>%
  arrange(year) %>%
  mutate(lagged_orangegonad = lag(mean_orangegonad)) 

scallopsummarydf_annualmeans <- na.omit(scallopsummarydf_annualmeans)

plot(scallopsummarydf_annualmeans$lagged_orangegonad, scallopsummarydf_annualmeans$mean_percentsmallrings)

# Make beta regression model and plot
mygam <- gam(mean_percentsmallrings~ lagged_orangegonad, family=betar(link="logit"), data = scallopsummarydf_annualmeans)
summary(mygam)
res1 <- simulateResiduals(mygam)
plot(res1)
min <- min(scallopsummarydf_annualmeans$lagged_orangegonad)
max <- max(scallopsummarydf_annualmeans$lagged_orangegonad)
new.x <- expand.grid(lagged_orangegonad = seq(min, max, length.out = 1000))
new.y <- predict(mygam, newdata = new.x, se.fit = TRUE, type="response")
new.y <- data.frame(new.y)
addThese <- data.frame(new.x, new.y)
addThese <- rename(addThese, y = fit, SE = se.fit)
addThese <- mutate(addThese, lwr = y - 1.96 * SE, upr = y + 1.96 * SE) 
addThese <- rename(addThese, mean_percentsmallrings = y)
RipePlot<-ggplot(scallopsummarydf_annualmeans, aes(x = lagged_orangegonad, y = mean_percentsmallrings)) +
  geom_point(size =2.5, alpha = .75)+
  geom_smooth(data = addThese, aes(ymin = lwr, ymax = upr), stat = 'identity',color="darkseagreen4")+
  theme_bw() +
  ylab("Percent Small Rings the Next Year")+
  xlab("Percent Ripe Scallops")+
  theme(text = element_text(size=10)) +
  theme(panel.background = element_blank())
RipePlot

ggsave("RipePlot.tiff",RipePlot, dpi = 300, bg = "white",
       width = 20,
       height = 15,
       units = "cm")

### Create scallopsummarydf_annualmeans with annual means for SEPERATE SITES

scallopsummarydf_annualmeans <- scallopringsummarydf %>%
  mutate(year = format(Date, "%Y")) %>%
  group_by(year, Site) %>%  
  summarise(
    mean_orangegonad = mean(orangegonad_percentripe, na.rm = TRUE),
    mean_percentsmallrings = mean(percentsmallrings_lessthan20mm, na.rm = TRUE)
  ) %>%
  ungroup()

scallopsummarydf_annualmeans <- scallopsummarydf_annualmeans %>%
  arrange(year) %>%
  mutate(lagged_orangegonad = lag(mean_orangegonad)) 

scallopsummarydf_annualmeans <- na.omit(scallopsummarydf_annualmeans)

plot(scallopsummarydf_annualmeans$lagged_orangegonad, scallopsummarydf_annualmeans$mean_percentsmallrings)

# Create the new 'location' column with grouped labels
scallopsummarydf_annualmeans<- scallopsummarydf_annualmeans %>%
  mutate(Site = case_when(
    Site %in% c("NW Harbor (Off Split Rock)", "NW Harbor - S of Alewife Creek", "NW Harbor - Barcelona Neck", "NW Harbor - N of Mile Hill Rd") ~ "NW Harbor",
    Site %in% c("Southold - off Cedar Beach", "Southold Bay - N central hole") ~ "Southold",
    Site %in% c("Robin's Island - W side", "Robin's Island - W side (dredged)") ~ "Robins",
    Site %in% c("Shelter Island - N Tip  (Off Hay Beach)", "Shelter Island - E side") ~ "Shelter Island",
    TRUE ~ Site 
  ))

mygam <- gam(mean_percentsmallrings~ lagged_orangegonad, family=betar(link="logit"), data = scallopsummarydf_annualmeans)
summary(mygam) 

#model has trouble converging when trying to keep sites/bays separate, so prob just go with previous model where we used annual means


####### OBJECTIVE: Examine relationship between the number of days past the start of fall and fecundity


scallopringsummarydf <- X_2024_master_scallop_ring_paper_summary_data_highlighting_revised15Mar2025_data_UNCHANGED_4

scallopringsummarydf$Date <- as.Date(scallopringsummarydf$Date)
scallopringsummarydf$dayspostsep30 <- as.numeric(scallopringsummarydf$Date - as.Date(paste0(format(scallopringsummarydf$Date, "%Y"), "-09-30")))

scallopringsummarydf <- scallopringsummarydf[!is.na(scallopringsummarydf$orangegonad_percentripe), ]
scallopringsummarydf$orangegonad_percentripe <- scallopringsummarydf$orangegonad_percentripe / 100

mygam <- gam(orangegonad_percentripe~ dayspostsep30, family=betar(link="logit"), data = scallopringsummarydf)
summary(mygam)
min <- min(scallopringsummarydf$dayspostsep30)
max <- max(scallopringsummarydf$dayspostsep30)
new.x <- expand.grid(dayspostsep30 = seq(min, max, length.out = 1000))
new.y <- predict(mygam, newdata = new.x, se.fit = TRUE, type="response")
new.y <- data.frame(new.y)
addThese <- data.frame(new.x, new.y)
addThese <- rename(addThese, y = fit, SE = se.fit)
addThese <- mutate(addThese, lwr = y - 1.96 * SE, upr = y + 1.96 * SE)
addThese <- rename(addThese, orangegonad_percentripe = y)
RipePlotDaysPost<-ggplot(scallopringsummarydf, aes(x = dayspostsep30, y = orangegonad_percentripe )) +
  geom_point(size =2.5, alpha = .75)+
  geom_smooth(data = addThese, aes(ymin = lwr, ymax = upr), stat = 'identity',color="darkseagreen4")+
  theme_bw() +
  ylab("Percent Ripe Scallops")+
  xlab("Days Post Sept 30")+
  theme(text = element_text(size=10)) +
  theme(panel.background = element_blank())
RipePlotDaysPost

ggsave("RipePlotDaysPost.tiff",RipePlotDaysPost, dpi = 300, bg = "white",
       width = 20,
       height = 15,
       units = "cm")


####### OBJECTIVE: Examine meat weight temporal trends


###### DENNIS CODE 


df<- read.csv("/Users/rayczaja/Desktop/Ring_Data_New_Dataset.csv")

# Convert date column to Date format
df$Date <- as.Date(df$Date, "%m/%d/%Y")

# Extract year from date
df$year <- format(df$Date, "%Y")

# Calculate number of days after September 30 for each year
df$days_post_sept30 <- ifelse(
  df$Date >= as.Date(paste(df$year, "-09-30", sep = "")),
  as.numeric(df$Date - as.Date(paste(df$year, "-09-30", sep = ""))),
  as.numeric(df$Date - as.Date(paste(as.numeric(df$year) - 1, "-09-30", sep = "")))
)

# Keep only samples between Sept 30 and Dec 31
df <- df[(format(df$Date, "%m-%d") >= "09-30" & format(df$Date, "%m-%d") <= "12-31"), ]

# Remove rows with missing values in key columns
df_clean <- df[complete.cases(df$days_post_sept30, df$meats_per_lbs, df$n), ]

# (Commented out) Option to require at least 2 data points per year
# df_clean <- df_clean %>%
#   group_by(year) %>%
#   filter(n() >= 2) %>%
#   ungroup()

# Convert year column to numeric and filter to years <= 2004
df_clean$year <- as.numeric(df_clean$year)
df_clean_pre <- df_clean %>% filter(year <= 2004)

# Filter out sites with fewer than 3 observations
df_clean_pre <- df_clean_pre %>%
  group_by(Site.Label) %>%
  filter(n() >= 3) %>%
  ungroup()

# Ensure year is numeric again (if modified above)
df_clean_pre$year <- as.numeric(df_clean_pre$year)

# Convert site labels to factor (categorical variable)
df_clean_pre$Site.Label <- as.factor(df_clean_pre$Site.Label)
levels(df_clean_pre$Site.Label)  # View levels (sites)

# Summarize number of observations per site
result_pre <- df_clean_pre %>%
  group_by(Site.Label) %>%
  summarise(`n=` = n(), .groups = 'drop')

# View site-wise sample sizes
result_pre

# Convert meats_per_lbs to weight per meat in grams and round
df_clean_pre$Meat_weight <- round((453.6 / df_clean_pre$meats_per_lbs), 1)

# Fit linear model: meat weight ~ days since Sept 30 + year + site
model_Day_pre <- lm(Meat_weight ~ days_post_sept30 + year + Site.Label, data = df_clean_pre)

# Create diagnostic plots for linear model (residuals, leverage, etc.)
par(mfrow = c(2, 2))  # 2x2 grid for model diagnostic plots
plot(model_Day, which = c(1:4))

# Perform ANOVA on model
ANV_pre <- anova(model_Day_pre)
ANV_pre  # Print ANOVA table

# Get full model summary with coefficients and significance
SUMM_pre <- summary(model_Day_pre)
SUMM_pre  # Print summary

# Plot meat weight by year, colored by site, with linear trend
plot1<-ggplot(df_clean_pre, aes(x = year, y = Meat_weight)) +
  geom_point(aes(size = n, color = Site.Label)) +
  geom_smooth(method = "lm", size =0.5, color ="black") +
  labs(
    y = "Average Meat Weight (g)",
    x = "Year",
    subtitle = "1990 through 2004"
  ) +
  theme_bw() +
  scale_size_continuous(name = "n", breaks = c(100, 300, 500)) +
  scale_color_manual(
    name = "Site",
    values = c("Hog Neck Bay" = "darkseagreen3", 
               "Shelter Island" = "indianred3", 
               "East Marion" = "skyblue3",
               "Orient Harbor" = "gold2",
               "Hallock Bay" = "plum3",
               "NW Harbor East Side" = "darkolivegreen4",
               "NW Harbor" = "slateblue2",
               "Hog Neck Bay" = "sienna3"))+
  guides(color = "none") +
  theme(
    legend.position = c(0.85, 0.80), 
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8))

# Plot meat weight by days after September 30th
plot2<-ggplot(df_clean_pre, aes(x = days_post_sept30, y = Meat_weight)) +
  geom_point(aes(size = n, color = Site.Label)) +
  geom_smooth(method = "lm", size =0.5, color ="black") +
  labs(
    y = "Average Meat Weight (g)",
    x = "Days Post September 30",
  ) +
  theme_bw() +
  scale_size_continuous(name = "n", breaks = c(100, 300, 500)) +
  scale_color_manual(
    name = "Site",
    values = c("Hog Neck Bay" = "darkseagreen3", 
               "Shelter Island" = "indianred3", 
               "East Marion" = "skyblue3",
               "Orient Harbor" = "gold2",
               "Hallock Bay" = "plum3",
               "NW Harbor East Side" = "darkolivegreen4",
               "NW Harbor" = "slateblue2",
               "Hog Neck Bay" = "sienna3"))+
  guides(color = "none") +
  theme(
    legend.position = c(0.85, 0.30), 
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8))

# Summarize sample counts again by site (for output)
result <- df_clean %>%
  group_by(Site.Label) %>%
  summarise(`n=` = n(), .groups = 'drop')

# Capture printed output of ANOVA, model summary, and site summary
output1 <- capture.output(print(ANV))
output2 <- capture.output(print(SUMM))
output3 <- capture.output(print(result))

# Define file paths to save outputs
file_path1 <- "C:/Users/Denni/OneDrive/Documents/Scallop_Ring_Project/Model_Day_Anova_1990to2004_04232025.txt"
file_path2 <- "C:/Users/Denni/OneDrive/Documents/Scallop_Ring_Project/Model_Day_Summary_1990to2004_04232025.txt"
file_path3 <- "C:/Users/Denni/OneDrive/Documents/Scallop_Ring_Project/n_per_site_1990to2004_04232025.txt"

# Write captured output to text files
writeLines(output1, con = file_path1)
writeLines(output2, con = file_path2)
writeLines(output3, con = file_path3)

############### 2005_through_2019 ########################################

# Load scallop ring dataset
df<- read.csv("/Users/rayczaja/Desktop/Ring_Data_New_Dataset.csv")

# Convert 'Date' column to Date format
df$Date <- as.Date(df$Date, "%m/%d/%Y")

# Extract year from the date and store as a character column
df$year <- format(df$Date, "%Y")

# Calculate number of days since September 30th of the given year
# If date is after Sept 30, difference from same-year Sept 30; else from prior year
df$days_post_sept30 <- ifelse(
  df$Date >= as.Date(paste(df$year, "-09-30", sep = "")),
  as.numeric(df$Date - as.Date(paste(df$year, "-09-30", sep = ""))),
  as.numeric(df$Date - as.Date(paste(as.numeric(df$year) - 1, "-09-30", sep = "")))
)

# Subset data to only include records between Sept 30 and Dec 31
df <- df[(format(df$Date, "%m-%d") >= "09-30" & format(df$Date, "%m-%d") <= "12-31"), ]

# Remove rows with missing values in key variables
df_clean <- df[complete.cases(df$days_post_sept30, df$meats_per_lbs, df$n), ]

# (Optional) You could filter out years with < 2 obs — commented out below
# df_clean <- df_clean %>%
#   group_by(year) %>%
#   filter(n() >= 2) %>%
#   ungroup()

# Convert 'year' column to numeric and filter to 2005 or later
df_clean$year <- as.numeric(df_clean$year)
df_clean_post <- df_clean %>% filter(year >= 2005)

# Only include sites with at least 3 observations
df_clean_post <- df_clean_post %>%
  group_by(Site.Label) %>%
  filter(n() >= 3) %>%
  ungroup()

# Reassert year as numeric (in case it was transformed)
df_clean_post$year <- as.numeric(df_clean_post$year)

# Convert 'Site.Label' column to factor (for modeling)
df_clean_post$Site.Label <- as.factor(df_clean_post$Site.Label)
levels(df_clean$Site.Label)  # View the included sites

# Summarize how many observations exist per site
result_post <- df_clean_post %>%
  group_by(Site.Label) %>%
  summarise(`n=` = n(), .groups = 'drop')
result_post  # Display summary

# Convert meats per pound to grams per individual meat, round to 1 decimal place
df_clean_post$Meat_weight <- round((453.6 / df_clean_post$meats_per_lbs), 1)

# Fit a linear model: meat weight as a function of time since Sept 30, year, and site
model_Day_post <- lm(Meat_weight ~ days_post_sept30 + year + Site.Label, data = df_clean_post)

# Show standard diagnostic plots for model fit
par(mfrow = c(2, 2))  # Arrange plots in a 2x2 grid
plot(model_Day, which = c(1:4))  # Residuals, leverage, etc.

# Perform ANOVA on the linear model
ANV_post<- anova(model_Day_post)
ANV  # Show ANOVA table

# Output full model summary (coefficients, p-values, R²)
SUMM_post <- summary(model_Day_post)
SUMM

# PLOT 4: Average meat weight by year
plot4<-ggplot(df_clean_post, aes(x = year, y = Meat_weight)) +
  geom_point(aes(size = n, color = Site.Label)) +
  geom_smooth(method = "lm", size =0.5, color ="black") +
  labs(
    y = "Average Meat Weight (g)",
    x = "Year",
    subtitle = "2005 through 2019"
  ) +
  theme_bw() +
  scale_size_continuous(name = "n", breaks = c(100, 200, 300)) +
  scale_color_manual(
    name = "Site",
    values = c("Hog Neck Bay" = "darkseagreen3", 
               "Noyack Bay" = "indianred3", 
               "Southold Bay" = "skyblue3"))+
  guides(color = "none") +
  theme(
    legend.position = c(0.85, 0.25), 
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    axis.title.y = element_blank(),
    axis.text.y = element_blank())

# PLOT 5: Average meat weight by days since Sept 30th
plot5<-ggplot(df_clean_post, aes(x = days_post_sept30, y = Meat_weight)) +
  geom_point(aes(size = n, color = Site.Label)) +
  geom_smooth(method = "lm", size =0.5, color ="black") +
  labs(
    y = "Average Meat Weight (g)",
    x = "Days Post September 30",
  ) +
  theme_bw() +
  scale_size_continuous(name = "n", breaks = c(100, 200, 300)) +
  scale_color_manual(
    name = "Site",
    values = c("Hog Neck Bay" = "darkseagreen3", 
               "Noyack Bay" = "indianred3", 
               "Southold Bay" = "skyblue3"))+
  guides(color = "none") +
  theme(
    legend.position = c(0.85, 0.25), 
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    axis.title.y = element_blank(),
    axis.text.y = element_blank())

DennisPlot<-grid.arrange(plot1, plot4, plot2, plot5,
                         ncol = 2, nrow = 2)

ggsave("DennisPlot.tiff",DennisPlot, dpi = 300, bg = "white",
       width = 20,
       height = 20,
       units = "cm")

# Recalculate sample size per site (for output logging)
result <- df_clean %>%
  group_by(Site.Label) %>%
  summarise(`n=` = n(), .groups = 'drop')

# Capture output text for logging
output1 <- capture.output(print(ANV))
output2 <- capture.output(print(SUMM))
output3 <- capture.output(print(result))

# Define file paths for saving analysis outputs
file_path1 <- "C:/Users/Denni/OneDrive/Documents/Scallop_Ring_Project/Model_Day_Anova_2005to2019_04232025.txt"
file_path2 <- "C:/Users/Denni/OneDrive/Documents/Scallop_Ring_Project/Model_Day_Summary_2005to2019_04232025.txt"
file_path3 <- "C:/Users/Denni/OneDrive/Documents/Scallop_Ring_Project/n_per_site_2005to2019_04232025.txt"

# Write output to text files
writeLines(output1, con = file_path1)
writeLines(output2, con = file_path2)
writeLines(output3, con = file_path3)



####### OBJECTIVE: Examine percent small ring relationship with temperature


temp_range <- range(annual_analyses$`temp_may-nov`, na.rm = TRUE)
pctnub_range <- range(annual_analyses$pctnub_yearsampled, na.rm = TRUE)

scale_factor <- diff(temp_range) / diff(pctnub_range)
offset <- temp_range[1] - pctnub_range[1] * scale_factor

ScallopTempPlot<-ggplot(annual_analyses, aes(x = year)) +
  annotate("rect", xmin = 2005, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "beige", alpha = 0.4) +
  geom_line(aes(y = `temp_may-nov`), color = "indianred3") +
  geom_point(aes(y = pctnub_yearsampled * scale_factor + offset),
             size = 2.5, alpha = 0.75) +
  scale_y_continuous(
    name = "Temperature (°C)",
    sec.axis = sec_axis(~ (. - offset) / scale_factor,
                        name = "Percent Small Rings"))+
  scale_x_continuous(name = "Year") +  
  annotate("text", x = 1993, y = max(annual_analyses$`temp_may-nov`, na.rm = TRUE),
           label = "Harvest Season Start: Oct 1", hjust = 0, vjust = -0.5, size = 4.2) +
  annotate("text", x = 2010, y = max(annual_analyses$`temp_may-nov`, na.rm = TRUE),
           label = "Harvest Season Start: Nov 1", hjust = 0, vjust = -0.5, size = 4.2) +
  geom_vline(xintercept = 2005, linetype = "dashed", color = "grey40") + 
  theme_bw()

ggsave("ScallopTempPlot.tiff",ScallopTempPlot, dpi = 300, bg = "white",
       width = 20,
       height = 15,
       units = "cm")



####### OBJECTIVE: Examine pea crab trends

#########
# this df takes the OG master df, but...
# a) removed the observation that had 'few' for the percent pea crab incidence 
# b) inserted 29 (the highest reported values) for the two observsations that had "lots" for the percent pea crab incidence
X2024_master_scallop_ring_paper_cleaned <- read_excel("Desktop/2024.master.scallop.ring.paper.cleaned.xls")

# remove all observations that have NA for percent incidence
peacrabdf <- X2024_master_scallop_ring_paper_cleaned %>%
  filter(!is.na(percent_incidence))

# conert to decimal for regression analyais
peacrabdf <- peacrabdf %>%
  mutate(percent_incidence = percent_incidence / 100)

# glms
m<-glm(percent_incidence ~ Year, binomial(link = "logit"), data=peacrabdf)
m<-glm.nb(percent_incidence ~ Year, data=peacrabdf)

# beta regression 
m<-betareg(percent_incidence ~ Year, data=peacrabdf)

summary(m)

# diagnostics for glms
res1 <- simulateResiduals(m)
plot(res1)

# diagnostics for beta reg
par(mfrow = c(3, 2))
suppressWarnings(RNGversion("3.5.0"))
set.seed(123)
plot(m, which = 1:4, type = "pearson")
plot(m, which = 5, type = "deviance", sub.caption = "")
plot(m, which = 1, type = "deviance", sub.caption = "")

# simple line plot
PeaCrabPlot<-ggplot(peacrabdf, aes(x = Year, y = percent_incidence)) +
  geom_point(size =2.5, alpha = .75)+
  geom_smooth(method = "lm", se = TRUE, color = "darkseagreen4") +
  labs(x = "Year",
       y = "Percent Incidence of Pea Crab (decimal conversion)") +
  theme_bw()

ggsave("PeaCrabPlot.tiff",PeaCrabPlot, dpi = 300, bg = "white",
       width = 20,
       height = 15,
       units = "cm")


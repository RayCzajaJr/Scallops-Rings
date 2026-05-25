
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
library(MASS)
library(glmmTMB)


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
    Location %in% c("NW Harbor (Barcelona Point)",  "NW Harbor (Barcelona Pt)", "NWHarbor (N of Mile Hill Road)") ~ "NW Harbor - Barcelona",
    Location == "Flanders" ~ "Flanders",
    Location == "Noyack Bay (EW side)" ~ "Noyack Bay - W Side",
    Location %in% c("Hallock Bay", "Hallock Bay (Outside Narrow River)", "Hallock Bay (central flats)") ~ "Hallock Bay",
    Location %in% c("Hog Neck", "Hog Neck Bay (SE corner)") ~ "Hog Neck",
    Location %in% c("OH Harbor", "OH Harbor (North)", "OH Harbor North", "OH Harbor - N") ~ "N Orient Harbor",
    Location %in% c("Robins Island (West central side)", "Robins Island (W side)") ~ "Robins Island",
    Location %in% c("Shelter Island (Hay Beach)") ~ "Shelter Is - Hay Beach",
    Location %in% c("Shelter Island (East Side)", "Shelter Island (NE Side)") ~ "Shelter Is - E Side",
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

sizedistscallops_nov <- subset(sizedistscallops_nov, Year != 1993)

# Setting up ridge plot aesthetics
year_levels <- sort(unique(c(sizedistscallops_oct$Year,
                             sizedistscallops_nov$Year)), decreasing = TRUE)

shared_xlim <- c(
  min(c(sizedistscallops_oct$`Total Ht`,
        sizedistscallops_nov$`Total Ht`), na.rm = TRUE),
  max(c(sizedistscallops_oct$`Total Ht`,
        sizedistscallops_nov$`Total Ht`), na.rm = TRUE))

nov_site_levels <- c(
  "Hallock Bay",
  "Southold - Cedar Beach",
  "Hog Neck",
  "Noyack Bay - W Side",
  "Robins Island",
  "Mattituck"
)

oct_site_levels <- c(
  "Hallock Bay",
  "N Orient Harbor",
  "Shelter Is - Hay Beach",
  "Shelter Is - E Side",
  "NW Harbor - E Side",
  "NW Harbor - Barcelona",
  "Hog Neck"
)

# Set up counts to be displayed on ridge plots
label_x <- shared_xlim[1] + diff(shared_xlim) * 0.02

oct_site_counts <- sizedistscallops_oct %>%
  filter(location %in% oct_site_levels) %>%
  group_by(location) %>%
  summarise(n_years = n_distinct(Year), .groups = "drop") %>%
  mutate(
    location_f = factor(location, levels = rev(oct_site_levels)),
    label = paste0(n_years, ifelse(n_years == 1, " year", " years")),
    x = label_x)

nov_site_counts <- sizedistscallops_nov %>%
  filter(location %in% nov_site_levels) %>%
  group_by(location) %>%
  summarise(n_years = n_distinct(Year), .groups = "drop") %>%
  mutate(
    location_f = factor(location, levels = rev(nov_site_levels)),
    label = paste0(n_years, ifelse(n_years == 1, " year", " years")),
    x = label_x)

oct_year_counts <- sizedistscallops_oct %>%
  group_by(Year) %>%
  summarise(n_sites = n_distinct(location), .groups = "drop") %>%
  mutate(
    Year_f = factor(Year, levels = year_levels),
    label = paste0(n_sites, ifelse(n_sites == 1, " site", " sites")),
    x = label_x)

nov_year_counts <- sizedistscallops_nov %>%
  group_by(Year) %>%
  summarise(n_sites = n_distinct(location), .groups = "drop") %>%
  mutate(
    Year_f = factor(Year, levels = year_levels),
    label = paste0(n_sites, ifelse(n_sites == 1, " site", " sites")),
    x = label_x)

# Ridge plot for shell and ring heights by site and year (two seperte months for shell heights)
novtotalhtyear <- ggplot(
  sizedistscallops_nov,
  aes(x = `Total Ht`,
    y = factor(Year, levels = year_levels),
    fill = factor(Year, levels = year_levels))) +
  stat_density_ridges( quantile_lines = TRUE,
    alpha = 0.7,
    quantiles = 2) +
  geom_text(data = nov_year_counts,
    aes(x = x, y = Year_f, label = label),
    inherit.aes = FALSE, hjust = 0, vjust= -0.5, size = 3) +
  scale_fill_carto_d(palette = "BurgYl") +
  scale_x_continuous(limits = shared_xlim) +
  theme_bw() +
  labs( x = NULL, y = NULL)+
  theme(legend.position = "none",
    plot.title = element_text(hjust = 0.5))

octtotalhtloc <- ggplot(
  sizedistscallops_oct,
  aes(x = `Total Ht`,
    y = factor(location, levels = rev(oct_site_levels)),
    fill = factor(location, levels = rev(oct_site_levels)))) +
  stat_density_ridges(quantile_lines = TRUE,
    alpha = 0.7, quantiles = 2) +
  geom_text(data = oct_site_counts,
    aes(x = x, y = location_f, label = label),
    inherit.aes = FALSE,
    hjust = 0,vjust= -0.5, size = 3) +
  scale_fill_carto_d(palette = "BluGrn") +
  scale_x_continuous(limits = shared_xlim) +
  theme_bw() +
  labs(x = NULL,y = NULL) +
  theme(legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_blank(),
    axis.text.y = element_text(angle = 45, hjust = 1))

octtotalhtyear <- ggplot(
  sizedistscallops_oct,
  aes( x = `Total Ht`,
    y = factor(Year, levels = year_levels),
    fill = factor(Year, levels = year_levels))) +
  stat_density_ridges(quantile_lines = TRUE,
    alpha = 0.7,
    quantiles = 2 ) +
  geom_text(data = oct_year_counts,
    aes(x = x, y = Year_f, label = label),
    inherit.aes = FALSE,hjust = 0,vjust= -0.5,size = 3) +
  scale_fill_carto_d(palette = "BurgYl") +
  scale_x_continuous(limits = shared_xlim) +
  theme_bw() +
  labs( x = NULL,y = NULL) +
  theme(legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_blank())

novtotalhtloc <- ggplot(
  sizedistscallops_nov,
  aes(x = `Total Ht`,
    y = factor(location, levels = rev(nov_site_levels)),
    fill = factor(location, levels = rev(nov_site_levels)))) +
  stat_density_ridges(quantile_lines = TRUE,
    alpha = 0.7,
    quantiles = 2) +
  geom_text(data = nov_site_counts,
    aes(x = x, y = location_f, label = label),
    inherit.aes = FALSE,hjust = 0,vjust= -0.5,size = 3) +
  scale_fill_carto_d(palette = "BluGrn") +
  scale_x_continuous(limits = shared_xlim) +
  theme_bw() +
  labs(x = NULL,y = NULL) +
  theme(legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_text(angle = 45, hjust = 1))

left_col <- arrangeGrob(
  octtotalhtloc,
  novtotalhtloc,
  ncol = 1,
  left = textGrob("Site", rot = 90, gp = gpar(fontsize = 16)))

right_col <- arrangeGrob(
  octtotalhtyear,
  novtotalhtyear,
  ncol = 1,
  left = textGrob("Year", rot = 90, gp = gpar(fontsize = 16)))

ridgeplots <- grid.arrange(
  left_col,
  right_col,
  ncol = 2,
  bottom = textGrob("Total Height", gp = gpar(fontsize = 16)))

ggsave("ridgeplots.tiff",ridgeplots, dpi = 300, bg = "white",
       width = 24,
       height = 24,
       units = "cm")


# Ring-height figure setup
combined_site_levels <- unique(c(oct_site_levels, nov_site_levels))
ring_xlim <- c(0, 80)
ring_label_x <- ring_xlim[2] - diff(ring_xlim) * 0.02

# Make a cleaned year column and include years from sizedistscallops too
sizedistscallops_ring <- sizedistscallops %>%
  mutate(Year_plot = as.numeric(as.character(Year)))

ring_year_levels <- sort(unique(c(year_levels, sizedistscallops_ring$Year_plot)), decreasing = TRUE)

ring_site_counts <- sizedistscallops_ring %>%
  filter(location %in% combined_site_levels) %>%
  group_by(location) %>%
  summarise(n_years = n_distinct(Year_plot, na.rm = TRUE), .groups = "drop") %>%
  mutate(location_f = factor(location, levels = rev(combined_site_levels)),
         label = paste0(n_years, ifelse(n_years == 1, " year", " years")),
         x = ring_label_x)

ring_year_counts <- sizedistscallops_ring %>%
  filter(!is.na(Year_plot)) %>%
  group_by(Year_plot) %>%
  summarise(n_sites = n_distinct(location), .groups = "drop") %>%
  mutate(Year_f = factor(Year_plot, levels = ring_year_levels),
         label = paste0(n_sites, ifelse(n_sites == 1, " site", " sites")),
         x = ring_label_x)

ringhtloc <- ggplot(
  sizedistscallops_ring %>% filter(location %in% combined_site_levels),
  aes(x = `Ring Ht`, y = factor(location, levels = rev(combined_site_levels)),
      fill = factor(location, levels = rev(combined_site_levels)))) +
  stat_density_ridges(quantile_lines = TRUE, alpha = 0.7, quantiles = 2) +
  geom_text(data = ring_site_counts, aes(x = x, y = location_f, label = label),
            inherit.aes = FALSE, hjust = 0.7, vjust = -0.5, size = 3) +
  scale_fill_carto_d(palette = "BluGrn") +
  scale_x_continuous(limits = ring_xlim) +
  theme_bw() +
  labs(x = NULL, y = NULL) +
  theme(legend.position = "none",
        axis.title = element_text(size = 16),
        axis.text.y = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))

ringhtyr <- ggplot(
  sizedistscallops_ring %>% filter(!is.na(Year_plot)),
  aes(x = `Ring Ht`, y = factor(Year_plot, levels = ring_year_levels),
      fill = factor(Year_plot, levels = ring_year_levels))) +
  stat_density_ridges(quantile_lines = TRUE, alpha = 0.7, quantiles = 2) +
  geom_text(data = ring_year_counts, aes(x = x, y = Year_f, label = label),
            inherit.aes = FALSE, hjust = 0.7, vjust = -0.5, size = 3) +
  scale_fill_carto_d(palette = "BurgYl") +
  scale_x_continuous(limits = ring_xlim) +
  theme_bw() +
  labs(x = NULL, y = NULL) +
  theme(legend.position = "none",
        axis.title = element_text(size = 16),
        plot.title = element_text(hjust = 0.5))

left_ring_col <- arrangeGrob(
  ringhtloc, left = textGrob("Site", rot = 90, gp = gpar(fontsize = 16)))

right_ring_col <- arrangeGrob(
  ringhtyr, left = textGrob("Year", rot = 90, gp = gpar(fontsize = 16)))

ringhts <- grid.arrange(
  left_ring_col, right_ring_col, ncol = 2, nrow = 1,
  bottom = textGrob("Ring Height", gp = gpar(fontsize = 16)))

ggsave("ringhts.tiff", ringhts, dpi = 300, bg = "white",
       width = 24, height = 20, units = "cm")

medians<-sizedistscallops %>%
  group_by(Year) %>%
  summarise(median_ring_ht = median(`Ring Ht`, na.rm = TRUE))


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

res <- residuals(m, type = "deviance")
acf(res, na.action = na.pass)  # Yup, significant temporal auto correlation. so let's make a new model that has a AR1 parameter

sizedistscallops_octandnov$Year_f <- factor(
  sizedistscallops_octandnov$Year,
  levels = sort(unique(sizedistscallops_octandnov$Year)))

m_ar1_year <- glmmTMB(Growthafterfirstwinter ~ Ring_Ht +
    (1 | location) + ar1(Year_f + 0 | group),
  data = transform(sizedistscallops_octandnov, group = factor(1)), REML = TRUE)

summary(m_ar1_year)
Anova(m)
r.squaredGLMM(m_ar1_year)

efct <- ggpredict(m_ar1_year, terms = "Ring_Ht [0:72 by=0.72]", type = "fixed")
efct <- as.data.frame(efct)

growthpostwinter <- ggplot() +
  geom_point(data = sizedistscallops_octandnov, aes(x = Ring_Ht, y = Growthafterfirstwinter),
    alpha = 0.4,color = "black") +
  geom_line(data = efct, aes(x = x, y = predicted),
    color = "darkseagreen4",linewidth = 1) +
  geom_ribbon(data = efct,aes(x = x, ymin = conf.low, ymax = conf.high),
    fill = "darkseagreen4",alpha = 0.3) +
  theme_bw() +
  labs(x = "Ring Height",y = "Growth after 1st Winter" ) +
  theme(plot.title = element_text(hjust = 0.5))
growthpostwinter

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
mygam <- gam(mean_percentsmallrings~ lagged_orangegonad, family=betar(link="logit"), data = scallopsummarydf_annualmeans) # annual means so no random effects 
summary(mygam)
res1 <- simulateResiduals(mygam)
plot(res1)
res <- residuals(mygam, type = "deviance")
acf(res, na.action = na.pass) # no evidence of temporal autocorrelation
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
  ylab("Proportion of Small Rings the Next Year")+
  xlab("Proportion of Ripe Scallops")+
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

scallopringsummarydf <- scallopringsummarydf %>%
  mutate(orangegonad_percentripe_beta = case_when(
      orangegonad_percentripe <= 0 ~ 0.001,
      orangegonad_percentripe >= 1 ~ 0.999,
      TRUE ~ orangegonad_percentripe))
mygam <- gam(orangegonad_percentripe_beta~ dayspostsep30, family=betar(link="logit"), data = scallopringsummarydf)
summary(mygam)
min <- min(scallopringsummarydf$dayspostsep30)
max <- max(scallopringsummarydf$dayspostsep30)
new.x <- expand.grid(dayspostsep30 = seq(min, max, length.out = 1000))
new.y <- predict(mygam, newdata = new.x, se.fit = TRUE, type="response")
new.y <- data.frame(new.y)
addThese <- data.frame(new.x, new.y)
addThese <- rename(addThese, y = fit, SE = se.fit)
addThese <- mutate(addThese, lwr = y - 1.96 * SE, upr = y + 1.96 * SE)
addThese <- rename(addThese, orangegonad_percentripe_beta = y)
RipePlotDaysPost<-ggplot(scallopringsummarydf, aes(x = dayspostsep30, y = orangegonad_percentripe_beta )) +
  geom_point(size =2.5, alpha = .75)+
  geom_smooth(data = addThese, aes(ymin = lwr, ymax = upr), stat = 'identity',color="darkseagreen4")+
  theme_bw() +
  ylab("Propoertion of Ripe Scallops")+
  xlab("Days Post Sept 30")+
  theme(text = element_text(size=10)) +
  theme(panel.background = element_blank())
RipePlotDaysPost

ggsave("RipePlotDaysPost.tiff",RipePlotDaysPost, dpi = 300, bg = "white",
       width = 20,
       height = 15,
       units = "cm")

# repeating the previous analysis but investigating site and year as random factors 

scallopringsummarydf <- scallopringsummarydf %>%
  mutate(Date = ymd(Date),
    Year = year(Date),
    Location = substr(Site, 1, 8))

scallopringsummarydf <- scallopringsummarydf %>%
  mutate(YearCat = as.factor(Year),Location = as.factor(Location))

unique(scallopringsummarydf$Location)

mygam <- gam(orangegonad_percentripe_beta~ dayspostsep30+
               s(Location, bs = "re") +
               s(YearCat, bs = "re"), family=betar(link="logit"), data = scallopringsummarydf)
summary(mygam)


#inspecting temporal auto correlation 
res <- residuals(mygam, type = "deviance")
acf(res, na.action = na.pass)

#this plot suggests no temporal autocorrelation, further further investigations (See Autocorrelation Investigation Script) revealed
# a) including an autoregressive term did not change the significance outcomes or smoothed shapes/effects
# b) the estimated AR parameter was small (0.1), and including this parameter decreased the quality of the fit

# repeating the predictions and plot to account for the random effects of year and location 
rng <- range(scallopringsummarydf$dayspostsep30, na.rm = TRUE)

new.x <- data.frame(
  dayspostsep30 = seq(rng[1], rng[2], length.out = 1000),
  Location = factor(levels(factor(scallopringsummarydf$Location))[1],
                    levels = levels(factor(scallopringsummarydf$Location))),
  YearCat = factor(levels(factor(scallopringsummarydf$YearCat))[1],
                   levels = levels(factor(scallopringsummarydf$YearCat)))
)

new.y <- predict(mygam, newdata = new.x, se.fit = TRUE, type = "response",
                 exclude = c("s(Location)", "s(YearCat)"))

addThese <- cbind(new.x, fit = new.y$fit, SE = new.y$se.fit) %>%
  mutate(lwr = fit - 1.96 * SE, upr = fit + 1.96 * SE)

RipePlotDaysPost <- ggplot(scallopringsummarydf,
                           aes(dayspostsep30, orangegonad_percentripe_beta)) +
  geom_point(size = 2.5, alpha = 0.75) +
  geom_ribbon(data = addThese,
              aes(x = dayspostsep30, ymin = lwr, ymax = upr),
              inherit.aes = FALSE, alpha = 0.25) +
  geom_line(data = addThese, aes(dayspostsep30, fit),
            color = "darkseagreen4", linewidth = 1) +
  theme_bw() +
  labs(x = "Days Post Sept 30", y = "Proportion of Ripe Scallops") +
  theme(text = element_text(size = 10),
        panel.background = element_blank())

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
model_Day_pre <- lmer(Meat_weight ~ days_post_sept30 + year + (1 | Site.Label),data = df_clean_pre)

# Create diagnostic plots for linear model (residuals, leverage, etc.)
par(mfrow = c(2, 2))  # 2x2 grid for model diagnostic plots
plot(model_Day_pre, which = c(1:4))

summary(model_Day_pre)

# Perform ANOVA on model
ANV_pre <- anova(model_Day_pre)
ANV_pre  # Print ANOVA table

# Get full model summary with coefficients, significance, etc
SUMM_pre <- summary(model_Day_pre)
SUMM_pre  
Anova(model_Day_pre, type = 3)
library(performance)
r2_nakagawa(model_Day_pre)
fixef(model_Day_pre)

# Check for temporal autocorrelation
res <- residuals(model_Day_pre, type = "deviance")
acf(res, na.action = na.pass)  

rename_sites <- c("East Marion" = "E Marion",
  "Hog Neck Bay" = "Hog Neck",
  "NW Harbor East Side" = "NW Harbor - E Side",
  "Orient Harbor" = "N Orient Harbor")

df_clean_pre$Site.Label <- dplyr::recode(df_clean_pre$Site.Label, !!!rename_sites)

# Plot meat weight by year, colored by site, with linear trend
pred_year <- ggpredict(
  model_Day_pre,
  terms = "year [all]",
  type = "fixed"
)

plot1 <- ggplot(df_clean_pre, aes(x = year, y = Meat_weight)) +
  geom_point(aes(color = Site.Label, size = n)) +
  geom_ribbon(
    data = pred_year,
    aes(x = x, ymin = conf.low, ymax = conf.high),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  geom_line(
    data = pred_year,
    aes(x = x, y = predicted),
    inherit.aes = FALSE,
    linewidth = 0.5,
    color = "black"
  ) +
  labs(
    y = "Average Meat Weight (g)",
    x = "Year",
    subtitle = "1990 through 2004"
  ) +
  theme_bw() +
  scale_size_continuous(name = "n", breaks = c(100, 300, 500)) +
  scale_color_manual(
    name = "Site",
    values = c(
      "Hog Neck" = "darkseagreen3", 
      "Shelter Island" = "indianred3", 
      "E Marion" = "skyblue3",
      "N Orient Harbor" = "gold2",
      "Hallock Bay" = "plum3",
      "NW Harbor - E Side" = "darkolivegreen4",
      "NW Harbor" = "slateblue2"),
    labels = c(
      "Hog Neck Bay" = "Hog Neck Bay",
      "Shelter Island" = "Shelter Island",
      "E Marion" = "E Marion",
      "N Orient Harbor" = "N Orient \nHarbor",
      "Hallock Bay" = "Hallock Bay",
      "NW Harbor - E Side" = "NW Harbor -\nE Side",
      "NW Harbor" = "NW Harbor",
      "Hog Neck" = "Hog Neck"
    )
  ) +
  theme(
    legend.position = "right", 
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )

plot1

pred_days <- ggpredict(
  model_Day_pre,
  terms = "days_post_sept30 [all]",
  type = "fixed"
)

# Plot meat weight by days after September 30th
plot2 <- ggplot(df_clean_pre, aes(x = days_post_sept30, y = Meat_weight)) +
  geom_point(aes(color = Site.Label, size = n)) +
  geom_ribbon(
    data = pred_days,
    aes(x = x, ymin = conf.low, ymax = conf.high),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  geom_line(
    data = pred_days,
    aes(x = x, y = predicted),
    inherit.aes = FALSE,
    linewidth = 0.5,
    color = "black"
  ) +
  labs(
    y = "Average Meat Weight (g)",
    x = "Days Post September 30"
  ) +
  theme_bw() +
  scale_size_continuous(name = "n", breaks = c(100, 300, 500)) +
  scale_color_manual(
    name = "Site",
    values = c(
      "Hog Neck" = "darkseagreen3", 
      "Shelter Island" = "indianred3", 
      "E Marion" = "skyblue3",
      "N Orient Harbor" = "gold2",
      "Hallock Bay" = "plum3",
      "NW Harbor - E Side" = "darkolivegreen4",
      "NW Harbor" = "slateblue2"
    ),
    labels = c(
      "Hog Neck Bay" = "Hog Neck Bay",
      "Shelter Island" = "Shelter Island",
      "E Marion" = "E Marion",
      "N Orient Harbor" = "N Orient \nHarbor",
      "Hallock Bay" = "Hallock Bay",
      "NW Harbor - E Side" = "NW Harbor -\nE Side",
      "NW Harbor" = "NW Harbor",
      "Hog Neck" = "Hog Neck"
    )
  ) +
  theme(
    legend.position = "right", 
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )

plot2

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
model_Day_post <- lmer(Meat_weight ~ days_post_sept30 + year + (1 | Site.Label),data = df_clean_post)

# Show standard diagnostic plots for model fit
par(mfrow = c(2, 2))  # Arrange plots in a 2x2 grid
plot(model_Day_post, which = c(1:4))  # Residuals, leverage, etc.

# Perform ANOVA on the linear model
ANV_post<- anova(model_Day_post)
ANV_post  # Show ANOVA table

# Output full model summary (coefficients, p-values, R²)
SUMM_post <- summary(model_Day_post)
SUMM_post
Anova(model_Day_post, type = 3)
r2_nakagawa(model_Day_post)
fixef(model_Day_post)


# Check for temporal autocorrelation
res <- residuals(model_Day_post, type = "deviance")
acf(res, na.action = na.pass)  # None, we're in the clear

rename_sites <- c(
  "Hog Neck Bay" = "Hog Neck"
)

df_clean_post$Site.Label <- dplyr::recode(df_clean_post$Site.Label, !!!rename_sites)

# PLOT 4: Average meat weight by year
pred_year_post <- ggpredict(model_Day_post,
  terms = "year [all]",
  type = "fixed")

pred_days_post <- ggpredict(model_Day_post,
  terms = "days_post_sept30 [all]",
  type = "fixed")

plot4 <- ggplot(df_clean_post, aes(x = year, y = Meat_weight)) +
  geom_point(aes(color = Site.Label, size = n)) +
  geom_ribbon(
    data = pred_year_post,
    aes(x = x, ymin = conf.low, ymax = conf.high),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  geom_line(
    data = pred_year_post,
    aes(x = x, y = predicted),
    inherit.aes = FALSE,
    linewidth = 0.5,
    color = "black"
  ) +
  labs(
    y = "Average Meat Weight (g)",
    x = "Year",
    subtitle = "2005 through 2019"
  ) +
  theme_bw() +
  scale_size_continuous(name = "n", breaks = c(100, 200, 300)) +
  scale_color_manual(
    name = "Site",
    values = c(
      "Hog Neck" = "darkseagreen3", 
      "Noyack Bay" = "hotpink3", 
      "Southold Bay" = "goldenrod3"
    )
  ) +
  theme(
    legend.position = "right", 
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    axis.title.y = element_blank(),
    axis.text.y = element_blank()
  )

plot4

plot5 <- ggplot(df_clean_post, aes(x = days_post_sept30, y = Meat_weight)) +
  geom_point(aes(color = Site.Label, size = n)) +
  geom_ribbon(
    data = pred_days_post,
    aes(x = x, ymin = conf.low, ymax = conf.high),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  geom_line(
    data = pred_days_post,
    aes(x = x, y = predicted),
    inherit.aes = FALSE,
    linewidth = 0.5,
    color = "black"
  ) +
  labs(
    y = "Average Meat Weight (g)",
    x = "Days Post September 30"
  ) +
  theme_bw() +
  scale_size_continuous(name = "n", breaks = c(100, 200, 300)) +
  scale_color_manual(
    name = "Site",
    values = c(
      "Hog Neck" = "darkseagreen3", 
      "Noyack Bay" = "hotpink3", 
      "Southold Bay" = "goldenrod3"
    )
  ) +
  theme(
    legend.position = "right", 
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    axis.title.y = element_blank(),
    axis.text.y = element_blank()
  )

plot5

DennisPlot<-grid.arrange(plot1, plot4, plot2, plot5,
                         ncol = 2, nrow = 2)

ggsave("DennisPlot.tiff",DennisPlot, dpi = 300, bg = "white",
       width = 28,
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
  annotate("text", x = 1992, y = max(annual_analyses$`temp_may-nov`, na.rm = TRUE),
           label = "Harvest Season: early October", hjust = 0, vjust = -0.5, size = 4.2) +
  annotate("text", x = 2009, y = max(annual_analyses$`temp_may-nov`, na.rm = TRUE),
           label = "Harvest Season: early November", hjust = 0, vjust = -0.5, size = 4.2) +
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
  labs(x = "Year",
       y = "Proportion of Scallops with Female Pea Crabs") +
  theme_bw()
PeaCrabPlot

ggsave("PeaCrabPlot.tiff",PeaCrabPlot, dpi = 300, bg = "white",
       width = 20,
       height = 15,
       units = "cm")


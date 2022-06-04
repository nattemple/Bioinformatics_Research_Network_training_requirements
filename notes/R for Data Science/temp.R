library(tidyverse)
library(ggplot2)
library(janitor)
library(readr)
library(tidyr)
library(plotly)
library(reshape2)
library(dplyr)
library(ggpubr)
library(kableExtra)

# 1.

my_data <- read_csv("gapminder_clean.csv")
my_data <- clean_names(my_data)

# 2.

s_columns <- select(
  my_data,
  "year",
  "gdp_percap",
  "co2_emissions_metric_tons_per_capita"
)
s_columns_1962_na <- filter(s_columns, year == 1962) %>% drop_na
s_columns_1962_na_no_outlier <- filter(s_columns, year == 1962 & gdp_percap < 75000) %>% drop_na
s_columns_na <- s_columns %>% drop_na


Plot1 <- ggplot(data = s_columns_1962_na) +
  geom_point(mapping = aes(
    x = co2_emissions_metric_tons_per_capita,
    y = gdp_percap
  ))

Plot1_no_outlier <- ggplot(data = s_columns_1962_na_no_outlier) +
  geom_point(mapping = aes(
    x = co2_emissions_metric_tons_per_capita,
    y = gdp_percap
  ))

# 3.

# check for normality of co2_emissions_metric_tons_per_capita and gdp_percap
ggplot(data = s_columns_1962_na, mapping = aes(x=gdp_percap)) + 
  geom_histogram(aes(y=..density..),fill="bisque",color="white",alpha=0.7) + 
  geom_density() +
  geom_rug() +
  labs(x='GDP per cap') +
  theme_minimal()

res3 <- cor.test(s_columns_1962_na$co2_emissions_metric_tons_per_capita,
                 s_columns_1962_na$gdp_percap,
                 method = "spearman"
)

# 4.

# corr by year
corr_by_year <- s_columns_na %>%
  group_by(year) %>%
  summarize(cor = cor(co2_emissions_metric_tons_per_capita, gdp_percap, method="spearman")) %>%
  arrange(desc(cor))

# q. 5
s5_columns <- select(
  my_data, "year",
  "gdp_percap",
  "co2_emissions_metric_tons_per_capita",
  "continent",
  "pop"
)
s5_columns_1967 <- filter(s5_columns, year == 1967)
s5_columns_1967_na <- drop_na(s5_columns_1967)


fig5 <- plot_ly(
  data = s5_columns_1967_na,
  x = ~co2_emissions_metric_tons_per_capita,
  y = ~gdp_percap,
  type = "scatter",
  color = ~continent,
  size = ~pop
) %>% layout(title = "CO2 emissions vs. GDP in 1967")

fig5_ <- ggplotly(fig5)

# q. #6

s6_columns <- select(
  my_data,
  "energy_use_kg_of_oil_equivalent_per_capita",
  "continent",
  "country_name",
  "pop",
  "year"
)

s6_columns <- drop_na(s6_columns)
s6_columns <- arrange(s6_columns, continent, year)

s6_summary <- s6_columns %>%
  group_by(continent, year) %>%
  summarise_at(
    vars(energy_use_kg_of_oil_equivalent_per_capita),
    list(energy_use_kg_of_oil_equivalent_per_capita = mean)
  )

Plot6 <- s6_summary %>%
  ggplot(aes(
    x = year,
    y = energy_use_kg_of_oil_equivalent_per_capita,
    group = continent, color = continent
  )) +
  geom_line()

s6_columns_1987 <- filter(s6_columns, year == 1987)

# Box Plot
Plot6_ <- ggboxplot(s6_columns_1987,
                    y = "energy_use_kg_of_oil_equivalent_per_capita",
                    x = "continent",
                    ylab = "Energy use (kg of oil equivalent per capita)",
                    xlab = FALSE
)


# q. #7

# Compute the analysis of variance
aov <- aov(energy_use_kg_of_oil_equivalent_per_capita ~ continent,
           data = s6_columns_1987
)
# Summary of the analysis
summary_aov <- summary(aov)

s7_columns <- select(
  my_data,
  "imports_of_goods_and_services_percent_of_gdp",
  "continent",
  "country_name",
  "year"
) %>%
  filter(year >= 1990) %>%
  drop_na()

s7_summary <- s7_columns %>%
  group_by(continent, year) %>%
  summarise_at(
    vars(imports_of_goods_and_services_percent_of_gdp),
    list(imports_of_goods_and_services_percent_of_gdp = mean)
  )


Plot7 <- s7_summary %>%
  ggplot(aes(
    x = year,
    y = imports_of_goods_and_services_percent_of_gdp,
    group = continent, color = continent
  )) +
  geom_line()

# Box Plot
s7_columns_europe_asia <- s7_columns %>%
  filter(continent == "Europe" | continent == "Asia")

Plot7_ <- ggboxplot(s7_columns_europe_asia,
                    y = "imports_of_goods_and_services_percent_of_gdp",
                    x = "continent",
                    ylab = "Imports of goods and services (% of GDP)",
                    xlab = "Continent"
)

# Compute the analysis of variance
aov7 <- aov(imports_of_goods_and_services_percent_of_gdp ~ continent,
            data = s7_columns_europe_asia
)
# Summary of the analysis
summary_aov7 <- summary(aov7)

# q. #8

s8_columns <- select(
  my_data,
  "population_density_people_per_sq_km_of_land_area",
  "continent",
  "country_name",
  "year"
) %>% drop_na()

s8_summary <-
  s8_columns %>%
  group_by(country_name) %>%
  summarise_at(
    vars(
      population_density_people_per_sq_km_of_land_area
    ),
    list(population_density_people_per_sq_km_of_land_area = mean)
  )


s8_summary_desc <- s8_summary %>%
  arrange(desc(population_density_people_per_sq_km_of_land_area))

# q. #9

# get first and last years for each country
s9_columns <- select(
  my_data,
  "life_expectancy_at_birth_total_years",
  "continent",
  "country_name",
  "year"
) %>% drop_na()

df <- s9_columns %>%
  select(country_name, year) %>%
  group_by(country_name) %>%
  filter(row_number() == 1 | row_number() == n())

df$first_last <- rep(c("first", "last"), times = nrow(df) / 2)
df <- df %>% pivot_wider(names_from = first_last, values_from = year)

# filter out those with no data before 1962
df <- filter(df, first <= 1962)

s9_columns <- subset(s9_columns, country_name %in% df$country_name)

# get earliest and latest life expectancy
df <- s9_columns %>%
  select(country_name, life_expectancy_at_birth_total_years) %>%
  group_by(country_name) %>%
  filter(row_number() == 1 | row_number() == n())

df$first_last <- rep(c("first", "last"), times = nrow(df) / 2)
df <- df %>% pivot_wider(
  names_from = first_last,
  values_from = life_expectancy_at_birth_total_years
)

df9 <- df %>%
  mutate(last = round(last)) %>%
  mutate(first = round(first)) %>%
  mutate(delta = last - first) %>%
  arrange(desc(delta))

names(df9)[names(df9) == "first"] <- "earliest"
names(df9)[names(df9) == "last"] <- "latest"

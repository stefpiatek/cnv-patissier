library(DBI)
library(RSQLite)
library(tidyverse)
library(knitr)
library(here)
library(broom)
library(purrr)
library(glue)
library(ggsci)
library(lubridate)
library(stringr)

capture <- "ICR_n2"

get_times <- function(capture){
  con <- dbConnect(SQLite(), dbname= here::here("../../output", glue("{capture}.sqlite")))
  
  callers <- con %>%
    tbl("callers")
  
  runs <- con %>%
    tbl("runs") %>%
    left_join(callers, by=c("caller_id" = "id"), suffix = c("", "_caller")) %>%
    collect() %>%
    mutate(set = capture)
  
  dbDisconnect(con)
  return(runs)
}

total_times <- map_df(.x = c("ICR", "ICR_n2", "ICR_n3"),
    .f = get_times) %>%
  mutate(duration = ymd_hms(duration) - ymd_hms("1970-01-01 00:00:00")) %>%
  mutate(name = str_remove(name, "_case"))

ggplot(total_times, aes(x = name, y = duration, colour = name)) +
  geom_boxplot() +
  labs(x = "CNV caller",
       y = "Run duration (minutes)") + 
  theme_bw() +
  scale_y_continuous(breaks=seq(0, 240, by = 30)) +
  ggsci::scale_color_jco() +
  guides(fill=FALSE, colour=FALSE) +
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  ggsave(here::here("../../output/plots_tables/duration.pdf"), height = 4, width = 4) +
  ggsave(here::here("../../output/plots_tables/duration.png"), height = 4, width = 4)

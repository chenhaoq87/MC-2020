library(readxl)
isolate_df  <- read_excel("~/Downloads/16S ST by Sample for Sarah 2020-02-21.xls")
micro_df <- read_excel("~/Downloads/Samantha Sample micro data 2020-02-21.xlsx")

# remove Sam R's rep post intervention data so only include VSL data
micro_df  <- micro_df %>%
  filter(`VSL Number` != "Aim2_Int1") %>%
  filter(`VSL Number` != "Aim2_Int2") %>%
  filter(`VSL Number` != "Aim2_Int3_Post") %>%
  filter(`VSL Number` != "Aim2_Int3_Pre") %>%
  filter(`VSL Number` != "Aim2_Prelim") 


# subset Upstate Rochester data
UR_isolate_df <- isolate_df %>%
  filter(Plant == "Upstate Rochester")
UR_micro_df <- micro_df %>%
  filter(Plant == "Upstate Rochester")


# merge df
# UR_df <- merge(UR_isolate_df,UR_micro_df,by=c("VSL Number","SampleID","Day"))

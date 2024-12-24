
install.packages("tidyr")
library(tidyr)

df <- data.frame(
  Sample1 = c("E. coli", "K. pneumonia"),
  Sample2 = c("n/a", "K. pneumoniae"),
  Gene = c("tetA", "blacmx")
)

df_long <- pivot_longer(df, cols = c(Sample1, Sample2), names_to = "Sample", values_to = "genus")
df_long <- df_long[df_long$genus != "n/a", ]

print(df_long)

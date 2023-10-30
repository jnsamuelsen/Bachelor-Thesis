## Tmp code ##
## Playground ##

# Create empty vector for storage of unique row names in expression matrix
unique_rownames <- c()

# Vector keeping track of rows which should be deleted
row_delete_tracker <- c()

# Counter for the for loop
counter <- 1

# Loading operators library
library(operators)

for (row in row.names(data_exprs)) {
  if (row %!in% unique_rownames) {
    append(x = unique_rownames, values = row)
  } else {
    append(x = row_delete_tracker, values = counter)
  }
  counter = counter+1
}


# Creating the data frame as a tibble
expression_data <- as_tibble(data_exprs)

# Adding gene symbols as a column to the tibble
add_column(expression_data, genename = 0)

# Remove duplicate row names # !!!
duplicated(rownames(data_exprs))

data_exprs[!duplicated(data_exprs$rownames),]



# Matrix scatter plot
ggplot <- function(...) ggplot2::ggplot(...) + scale_color_brewer(palette = "Spectral")

ggpairs(tibble_data, columns = c("T.cells.CD8", "T.cells.CD4.naive")) +
  ggplot(data = tibble_data, mapping = aes(color = status))



# Finding the indices in the expression data frame where the uniquely chosen largest probe ids are 
largest_probeids <- (which(row.names(df_exprs)%in%unique_ids))
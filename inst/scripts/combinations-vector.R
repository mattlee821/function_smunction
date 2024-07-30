# dataframe with all (independent) combinations of vectors
a <- c("abc", "cde", "def")

# Calling expand.grid() Function
a <- expand.grid(x1, x1)

# remove x on x
a <- a[a$Var1 != a$Var2,]

# remove duplicates (reversed)
fltr <- !duplicated(apply(a, 1, function(x) paste0(sort(x), collapse = "")))
a <- a[fltr, ]

println <- function(...) {

  cat(date(), '>', ..., '\n')
  flush.console()

}

get_cores <- function(...) ifelse(interactive(), 1, future::availableCores())

mclapply2 <- function(X, FUN, ...,
                      mc.preschedule = TRUE, mc.set.seed = TRUE,
                      mc.silent = FALSE, mc.cores = get_cores(),
                      mc.cleanup = TRUE, mc.allow.recursive = FALSE) {

  original_class <- class(X)

  if ( is.data.frame(X) ) new_x <- FUN(X, ...)
  else                    new_x <- parallel::mclapply(X, FUN, ...,
                                                      mc.preschedule = mc.preschedule,
                                                      mc.set.seed = mc.set.seed,
                                                      mc.silent = mc.silent, mc.cores = mc.cores,
                                                      mc.cleanup = mc.cleanup,
                                                      mc.allow.recursive = mc.allow.recursive)

  class(new_x) <- original_class

  return ( new_x )

}

get_file_inputs <- function(command_args, ...) {

  inputs <- list(data=command_args[6], annotation=command_args[7])
  stopifnot(! is.na(inputs$data) | ! file.exists(inputs$data) | ! is.na(inputs$data))

  return ( inputs )

}

get_controls <- function(command_args, ...) {

  controls <- list(positive=command_args[9], negative=command_args[8])
  if ( all(!is.na(controls)) ) println('Controls are', controls$positive, '(positive) and',
                                       controls$negative, '(negative)')

  return ( controls )

}

check_headers <- function(dataframe,
                          essential_headers=c('compound', 'concentration', 'strain', 'plate_name', 'count'))
  stopifnot(all(essential_headers %in% names(dataframe)))

pyjoin <- function(x, char) paste(x, collapse=char)

splitter <- function(x, char, n=NULL, rejoin_char=char) {

  vapply(x, function(y) {

    split_string <- strsplit(y, char, fixed=TRUE)[[1]]

    if ( is.null(n) ) n_ <- seq_len(length(split_string)) else n_ <- n

    pyjoin(split_string[n_], char=rejoin_char)

  }, 'a')

}

multigrepl <- function(patterns, targets) apply(vapply(patterns, grepl, grepl(patterns[1], targets), targets), 1, any)

capitalize <- function(x) {

  initials <- substr(x, 1, 1)
  other_chars <- substr(x, ifelse(nchar(x) > 1, 2, 0), ifelse(nchar(x) > 1, nchar(x), 0))

  return ( paste0(toupper(initials), other_chars) )

}

#' @importFrom magrittr %>%
column2vector <- function (x, column) x %>% dplyr::ungroup() %>% dplyr::select_(column) %>% unlist()

#' @importFrom magrittr %>%
get_unique_values <- function(x, column, output.class=NULL) {

  x2 <- column2vector(x, column)

  if (inherits(x2, 'factor') ) x2 <- x2 %>% droplevels() %>% levels()
  else                         x2 <- unique(x2)

  if ( !is.numeric(x2) & !is.null(output.class) & is.function(output.class) )  x2 <- output.class(x2)

  return ( x2 )

}


#' @importFrom magrittr %>%
get_compound_core <- function(x, column='compound') {

  library('dplyr')

  new_df <- x %>%
    dplyr::ungroup() %>%
    dplyr::select_(column) %>%
    dplyr::distinct() %>%
    dplyr::arrange() %>%
    dplyr::mutate_(compound_stem=paste0('sapply(as.character(', column, '), function(x) strsplit(x, "-", fixed=TRUE)[[1]][2])'))

  return ( inner_join(x, new_df) )

}

#' @importFrom magrittr %>%
make_columns_factors <- function(x, columns) {

  mutate_string <- paste0('as.factor(', columns, ')') %>% setNames(columns)

  x <- x %>% dplyr::mutate_(.dots=mutate_string)

  return ( x )

}

check_same_levels <- function(df1, df2, column) {

  levels_present1 <- get_unique_values(df1, column)
  levels_present2 <- get_unique_values(df2, column)

  return(length(setdiff(levels_present1, levels_present2)) > 0)

}

# statistical functions
optimizer <- function(f, interval, plot=NULL,  ..., lower = min(interval),
                      upper = max(interval), maximum = FALSE,
                      tol = .Machine$double.eps^0.25, max_recursion=100, recursion_no=0) {

  optimized <- optimize(f=f, interval=interval, ..., lower=lower,
                        upper=upper, maximum=maximum, tol=tol)

  if ( is.finite(optimized$objective) & f(upper + 0.1) > optimized$objective | recursion_no > max_recursion ) {

    return ( optimized )

  } else {

    println('Shifting window from',
            pyjoin(signif(interval, 2), char=', '),
            'to',
            pyjoin(signif(interval, 2) + 1, char=', '),
            '; current minimum is',
            signif(optimized$objective, 2), 'at',
            signif(optimized$minimum, 2))
    optimizer(f=f, interval=interval + 1, plot=plot, ..., maximum=maximum, tol=tol,
              max_recursion=max_recursion, recursion_no=recursion_no + 1)

  }

}


cr_adjustment <- function(means, dispersion, model_matrix) {

  weights_matrix <- diag(1 / ((1 / means) + dispersion))
  fisher_info <- crossprod(model_matrix, weights_matrix) %*% model_matrix

  cr_penalty <- 0.5 * determinant(fisher_info, logarithm=TRUE)$modulus[1]

  return(cr_penalty)

}

nb_dispersion_log_likelihood <- function(observed_data, means, dispersion)
  sum(log(mapply(dnbinom, x=observed_data, mu=means, size=1 / dispersion)))


cr_adjusted_log_likelihood <- function(observed_data, means, dispersion, model_matrix)
  nb_dispersion_log_likelihood(observed_data,
                               means,
                               dispersion) - cr_adjustment(means,
                                                           dispersion,
                                                           model_matrix)

#' @importFrom errR %except%
estimate_nb_dispersion_mle <- function(model, data, initial_guess,
                                    log_likelihood_function=cr_adjusted_log_likelihood,
                                    bounds=c(-3, 3)) {

  model_string      <- as.character(model)

  model_glm <- glm(model, MASS::negative.binomial(1 / initial_guess), data) %except% NULL

  if ( is.null(model_glm) ) {

    println('Intial dispersion model didn\'t fit:', model_string)
    print(initial_guess)

  }

  response_variable <- names(model_glm)[1]
  response_vector   <- model_glm$model[ , 1]

  predicted_mus <- model_glm$fitted.values
  model_matrix  <- model.matrix(model_glm)

  objective_f <- function(x)
    -log_likelihood_function(observed_data=response_vector,
                             means=predicted_mus, dispersion=x, model_matrix=model_matrix)

  println('Optimizing log likelihood function')
  #optimized            <- optimizer(objective_f, bounds + log(initial_guess))
  optimized            <- nlm(objective_f, initial_guess)
  #print(optimized)
  optimized_dispersion <- optimized$estimate

  return ( optimized_dispersion )

}

hijack <- function (f, ...) {
  .f <- f
  args <- list(...)
  .f.formals <- formals(.f)
  new.formals <- names(.f.formals) %>% lapply(function(x) {
    if (x %in% names(args)) return(args[[x]]) else return(.f.formals[[x]])
  }) %>% setNames(names(.f.formals))
  #print(new.formals)
  formals(.f) <- new.formals
  return ( .f )
}


trim_glm <- function(x, ...) {

  x$data <- NULL
  x$y <- NULL
  #x$linear.predictors <- NULL
  x$weights <- NULL
  x$fitted.values <- NULL
  x$model <- NULL
  x$prior.weights <- NULL
  x$residuals <- NULL
  x$effects <- NULL
  #x$qr$qr <- NULL

  return ( x )

}

#' @importFrom magrittr %>%
calculate_roc_stats <- function(x) {

  x <- x %>%
    dplyr::summarize(n_positive=sum(condition_positive, na.rm=TRUE),
                     n_negative=sum(condition_negative, na.rm=TRUE),
                     n_tests=n(),

                     prevalence=mean(condition_positive, na.rm=TRUE),

                     n_called_positive=sum(called_positive, na.rm=TRUE),
                     n_called_negative=sum(called_negative, na.rm=TRUE),

                     n_false_positive=sum(false_positive, na.rm=TRUE),
                     n_false_negative=sum(false_negative, na.rm=TRUE),
                     n_true_positive=sum(true_positive, na.rm=TRUE),
                     n_true_negative=sum(true_negative, na.rm=TRUE),

                     true_positive_rate=sum(true_positive, na.rm=TRUE) / n_positive,
                     false_positive_rate=sum(false_positive, na.rm=TRUE) / n_negative,
                     true_negative_rate=sum(true_negative, na.rm=TRUE) / n_negative,
                     false_negative_rate=sum(false_negative, na.rm=TRUE) / n_positive,

                     accuracy=mean(correct, na.rm=TRUE),
                     false_discovery_rate=sum(false_positive, na.rm=TRUE) / n_called_positive,
                     positive_predictive_value=sum(true_positive, na.rm=TRUE) / n_called_positive,
                     f1_score=2 / ((1 / true_positive_rate) + (1 / positive_predictive_value))) %>%
    mutate(false_discovery_rate=ifelse(is.na(false_discovery_rate), 0, false_discovery_rate),
           positive_predictive_value=ifelse(is.na(positive_predictive_value), 0, positive_predictive_value),
           f1_score=ifelse(is.na(f1_score), 0, f1_score))

  return ( x )
}

#' @importFrom magrittr %>%
calculate_roc_table <- function(x, positive_compound, cutoff_column, cutoffs, prevalence) {

  n_positives <- nrow(dplyr::filter(x, grepl(positive_compound, compound)))

  n_negatives <- nrow(dplyr::filter(x, !grepl(positive_compound, compound)))

  n_to_sample <- prevalence * n_negatives

  sampled_positives <- x %>%
    dplyr::filter(grepl(positive_compound, compound)) %>%
    dplyr::ungroup() %>%
    dplyr::sample_n(n_to_sample)

  x2 <- x %>%
    dplyr::filter(!grepl(positive_compound, compound))
    dplyr::bind_rows(sampled_positives)

  x <- x %>%
    tidyr::crossing(cutoff=cutoffs) %>%
    dplyr::mutate(condition_positive=grepl(positive_compound, compound),
                  condition_negative=! condition_positive) %>%
    dplyr::mutate_(.dots=c(called_positive=paste(cutoff_column, '<= cutoff'))) %>%
    dplyr::mutate(called_negative=!called_positive,
                  true_positive=condition_positive & called_positive,
                  false_positive=condition_negative & called_positive,
                  true_negative=condition_negative & called_negative,
                  false_negative=condition_positive & called_negative,
                  correct=true_positive | true_negative)

  return ( x )

}



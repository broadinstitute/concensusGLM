#' @title Get rough Negative Binomial dispersion parameters
#' @description Get rough Negative Binomial dispersion parameters using ordinary least squares.
#' @details Based on the assumption:
#' \deqn{Var(K) = \mu + \alpha  \mu^2}
#' where K is count, \eqn{\mu} is the mean count, and \eqn{\alpha} is the Negative Binomial dispersion parameter. This is rearranged to:
#' \deqn{\alpha = (Var(K) - \mu) / \mu^2}
#' The result is usually way off the final dispersion, but at least gives something reasonable to work from when estimating batch
#' effects using a Negative Binomial GLM.
#' @param x concensusWorkflow or concensusDataSet.
#' @param ... Other arguments.
#' @return concensusWorkflow or concensusDataSet with a new \code{mean_variance_relationship} and \code{dispersion} attribute.
#' @export
getRoughDispersions <- function(x, ...) UseMethod('getRoughDispersions')

#' @rdname getRoughDispersions
#' @export
getRoughDispersions.default <- function(x, ...) stop('Can\'t get rough dispersions of ', class(x))

#' @rdname getRoughDispersions
#' @export
getRoughDispersions.concensusWorkflow <- function(x, ...) {

  x <- workflows::delay(x, getRoughDispersions, ...)

  return (x)

}

#' @rdname getRoughDispersions
#' @importFrom magrittr %>%
#' @export
getRoughDispersions.concensusDataSet <- function(x, grouping=c('compound', 'concentration', 'strain'), ...) {

  println('Calculating rough dispersions')

  x$mean_variance_relationship <- x$data %>%
    dplyr::group_by_(.dots=grouping) %>%
    dplyr::summarize(var_count=var(count),
              mean_count=mean(count)) %>%
    dplyr:: ungroup() %>%
    dplyr::mutate(mean_count_sq=mean_count^2)

  # NB model: var = mu + alpha * mu^2

  dispersion_model <- lm(var_count ~ 0 + mean_count_sq + offset(mean_count), x$mean_variance_relationship)

  dispersion <- coef(dispersion_model)[['mean_count_sq']]

  x$dispersion <- data.frame(rough_dispersion=dispersion)

  return ( x )

}

#' @title Get experimental batch effects
#' @description Estimate experimental batch effects.
#' @details Columns in the \code{data} attribute of the \code{concensusDataSet} object which
#' have fewer than 100 unique values and which are not named according to a list of stopwords are assumed to be experimental
#' handling annotations, and will be modeled using a Negative Binomial GLM, using dispersion parameters in the \code{dispersion}
#' attribute.
#'
#' Only the negative control observations will be used for estimating batch effects.
#'
#' A new column, \code{predicted_null_count}, is added to \code{data} of \code{concensusDataSet} with a prediction of
#' number of counts for every observation, assuming only batch effects (and no real effect of interest).
#'
#' If the model fitting fails, rather than throw off the whole pipeline, a \code{predicted_null_count} value of \code{1} is used.
#' @param x concensusWorkflow or concensusDataSet.
#' @param ... Other arguments.
#' @return concensusWorkflow or concensusDataSet with a new \code{batch_effect_model} and a new \code{predicted_null_count} column
#' in the \code{data} attribute.
#' @seealso \link{glm}
#' @export
getBatchEffects <- function(x, ...) UseMethod('getBatchEffects')

#' @rdname getBatchEffects
#' @export
getBatchEffects.default <- function(x, ...) stop('Can\'t get batch effects of ', class(x))

#' @rdname getBatchEffects
#' @export
getBatchEffects.concensusWorkflow <- function(x, ...) {

  x <- workflows::delay(x, getBatchEffects, ...)

  return (x)

}

#' @rdname getBatchEffects
#' @importFrom magrittr %>%
#' @importFrom errR %except%
#' @export
getBatchEffects.concensusDataSet <- function(x, grouping=c('compound', 'concentration', 'strain'), ...) {

  not_batch_effects <- unique(c('plate_number', 'column', 'plate_seq', 'compound', 'concentration',
                         'strain', 'solvent', 'replicate', 'plate_seq', 'plate_source',
                         'positive_control', 'negative_control', 'predicted_null_count', 'count',
                         'plate_type', 'well_type', 'well', 'quadrant', 'gc_content', 'condition_group',
                         'i5_id', 'i7_id'), grouping)

  untreated_data <- x$data %>% dplyr::filter(negative_control)

  println('Using', nrow(untreated_data), 'negative control observations to calculate batch effects')

  candidate_batch_effects <- names(x$data)[!names(x$data) %in% not_batch_effects]

  #print(candidate_batch_effects)

  selected_batch_effects <- Filter(function(z) {

    levels_present_untreated <- get_unique_values(untreated_data, z)
    #print(c(z, levels_present_untreated, ))
    return(length(levels_present_untreated) > 1 & length(levels_present_untreated) < 100)
  },
  candidate_batch_effects)

  println('Modeling batch effects:', pyjoin(selected_batch_effects, ', '))

  x$data <- x$data %>% make_columns_factors(selected_batch_effects)
  untreated_data <- x$data %>% dplyr::filter(negative_control)

  underrepresented_batch_effects <- Filter(function(z) check_same_levels(x$data, untreated_data, z), selected_batch_effects)

  if ( length(underrepresented_batch_effects) > 0 )
    println('Batch effects', pyjoin(underrepresented_batch_effects, ', '), 'are underrepresented; resampling')

  for ( be in underrepresented_batch_effects ) {

    levels_present_untreated <- get_unique_values(untreated_data, be, as.character)
    levels_present_treated   <- get_unique_values(x$data, be, as.character)

    missing_levels <- setdiff(levels_present_treated, levels_present_untreated)

    not_represented_yet    <- length(missing_levels) > 0

    average_representation <- untreated_data %>%
      dplyr::group_by_(.dots=be) %>%
      dplyr::summarize(n_times=n()) %>%
      dplyr::select(n_times) %>%
      unlist() %>% mean() %>% as.integer()

    if ( not_represented_yet ) {

      println(be, pyjoin(missing_levels, ', '),
              'with no untreated wells; sampling', average_representation,
              'wells as untreated')

      untreated_data <- bind_rows(untreated_data,
                                  x$data %>%
                                    dplyr::filter(!negative_control, !positive_control) %>%
                                    dplyr::filter_(paste0(be, '%in% c(\'', pyjoin(missing_levels, '\', \''), '\')')) %>%
                                    dplyr::sample_n(average_representation, replace=TRUE))
    }

  }

  underrepresented_batch_effects <- Filter(function(z) check_same_levels(x$data, untreated_data, z), selected_batch_effects)
  if ( length(underrepresented_batch_effects) > 0 )
    stop('Missing batch effect levels: ', pyjoin(as.character(underrepresented_batch_effects), ', '))

  model_formula <- as.formula(paste('count ~', pyjoin(selected_batch_effects, '+')))

  x$batch_effect_model <- glm(model_formula,
                              MASS::negative.binomial(1 / x$dispersion$rough_dispersion), untreated_data) %except% NA

  if ( !is.na(x$batch_effect_model) &  !is.null(x$batch_effect_model) )
    x$data <- x$data %>% dplyr::mutate(predicted_null_count=predict(x$batch_effect_model, ., type='response'))
  else
    x$data <- x$data %>% dplyr::mutate(predicted_null_count=1)

  return ( x )

}

#' @title Estimate Negative Binomial dispersion paramteter taking into account batch effects
#' @description Estimate Negative Binomial dispersion paramteter taking into account experimental batch effects.
#' @details Uses the CR penalized maximum profile likelihood method, holding the \eqn{\mu} of a GLM fixed and finding the
#' optimal dispersion \eqn{\alpha} using a Newton-type algorithm as implemented in \code{nlm}.
#'
#' If the \code{predicted_null_count} column is present in the \code{data} attribute of \code{concensusDataSet}, it is added to
#' the GLM as an \code{offset}. If the batch effects are real, this should raise the likelihood of and shrink the size of
#' the final dispersion parameter.
#'
#' This method will find a dispersion value with and without taking into account \code{predicted_null_count}, saving both results to
#' the columns of the \code{dispersion} attribute of \code{concensusDataSet}.
#' @param x concensusWorkflow or concensusDataSet.
#' @param ... Other arguments.
#' @return concensusWorkflow or concensusDataSet with a new \code{small_model_dispersion} and a new
#' \code{full_model_dispersion} column in the \code{dispersion} attribute.
#' @seealso \link{nlm}, \link{glm}
#' @export
getFinalDispersions <- function(x, ...) UseMethod('getFinalDispersions')

#' @rdname getFinalDispersions
#' @export
getFinalDispersions.default <- function(x, ...) stop('Can\'t get final dispersions of ', class(x))

#' @rdname getFinalDispersions
#' @export
getFinalDispersions.concensusWorkflow <- function(x, ...) {

  x <- workflows::delay(x, getFinalDispersions, ...)

  return (x)

}

#' @rdname getFinalDispersions
#' @importFrom magrittr %>%
#' @export
getFinalDispersions.concensusDataSet <- function(x, ...) {

  println('Getting final dispersions...')

  small_model_formula <- count ~ 1
  full_model_formula  <- count ~ 0 + log(predicted_null_count)

  #print(head(x$data))

  untreated_data         <- x$data %>% dplyr::filter(negative_control)

  small_model_dispersion_ <- estimate_nb_dispersion_mle(initial_guess=x$dispersion$rough_dispersion,
                                                        data=untreated_data,
                                                        model=small_model_formula)

  x$dispersion <- x$dispersion %>%
    dplyr::mutate(small_model_dispersion=small_model_dispersion_)

  if ( 'predicted_null_count' %in% names(untreated_data)) {

    full_model_dispersion_  <- estimate_nb_dispersion_mle(initial_guess=x$dispersion$rough_dispersion,
                                                          data=untreated_data,
                                                          model=full_model_formula)

    #print()

    x$dispersion <- x$dispersion %>%
      dplyr::mutate(full_model_dispersion=full_model_dispersion_)

  }

  #print(x$dispersion)

  return ( x )

}

#' @title Fit final model
#' @description Fit the final Negative Binomial GLM, taking into account batch effects and final dispersion estimates.
#' @details Using \code{predicted_null_count} column in the \code{data} attribute of \code{concensusDataSet} is added to
#' the GLM as an \code{offset}. A new column, \code{condition_group} is created as a concatenation of the combinations of values
#' encountered in the columns specified in \code{conditions}, which is used as a catagorical variable.
#'
#' The \code{condition_group} associated with the negative control is identified, and set as the reference. A Negative Binomial
#' GLM is fitted using no intercept, \code{predicted_null_count} as an \code{offset}, \code{condition_group} as a predictor.
#' The \code{count} column is the response variable.
#'
#' Since a \code{log} link is used, and \code{predicted_null_count} is estimated from native control batch effects, the extracted
#' coefficient estimates can be intreprested as log(fold change) realtive to the negative control.
#' @param x concensusWorkflow or concensusDataSet.
#' @param conditions Character vector. Columns in the \code{data} attribute of \code{concensusDataSet} which together identify groups
#' of replicates of the same condition.
#' @param grouping Character vector. Columns in the \code{data} attribute of \code{concensusDataSet} which together identify
#' analytically independent chunks of data, e.g. strains.
#' @param ... Other arguments.
#' @return concensusWorkflow or concensusDataSet with a new \code{model_parameters} attribute containing effect sizes (LFC) and
#' p-values.
#' @seealso \link{glm}
#' @export
getFinalModel <- function(x, ...) UseMethod('getFinalModel')

#' @rdname getFinalModel
#' @export
getFinalModel.default <- function(x, ...) stop('Can\'t get final model of ', class(x))

#' @rdname getFinalModel
#' @export
getFinalModel.concensusWorkflow <- function(x, ...) {

  x <- workflows::delay(x, getFinalModel, ...)

  return (x)

}

#' @rdname getFinalModel
#' @importFrom magrittr %>%
#' @import methods
#' @export
getFinalModel.concensusDataSet <- function(x, conditions=c('compound', 'concentration'), grouping='strain',
                                           ...) {

  println('Fitting final model...')

  x$data <- x$data %>%
    dplyr::ungroup() %>%
    dplyr::mutate_(.dots=c(condition_group=paste0('paste(', pyjoin(conditions, ', '), ', sep="__")')))

  negative_controls <- x$data %>% dplyr::filter(negative_control)

  negative_level <- (negative_controls %>% get_unique_values('condition_group', sort))[1]

  nb_dispersion <- x$dispersion$full_model_dispersion[1]

  println('Reference level is', negative_level)
  println('Dispersion is', signif(nb_dispersion, 2))
  println('Fitting Negative Binomial GLM for',
          x$data %>% dplyr::filter(!negative_control) %>% get_unique_values('condition_group', length),
          'conditions as defined', pyjoin(conditions, ' + '))

  x$model_parameters <- x$data %>%
    dplyr::filter(!negative_control) %>%
    dplyr::group_by_(.dots=c(grouping, conditions, 'condition_group')) %>%
    dplyr::do((glm(count ~ 0 + offset(log(predicted_null_count)) + condition_group,
                   family=MASS::negative.binomial(1 / nb_dispersion),
                   data= dplyr::bind_rows(., negative_controls) %>%
                     dplyr::mutate(condition_group=factor(condition_group) %>% relevel(ref=negative_level)),
                   y=FALSE) %>%
                 broom::tidy() %>%
                 dplyr::mutate(p.value=as.numeric(p.value),
                               lfc=estimate,
                               l2fc=estimate / log(2),
                               std.error_log2=std.error / log(2))) ) #%except% data.frame(NULL))

  intercepts <- negative_controls %>%
    dplyr::group_by_(.dots=c(grouping)) %>%
    dplyr::summarize(mean_intercept=mean(predicted_null_count, na.rm=TRUE),
                     var_intercept=var(predicted_null_count, na.rm=TRUE))

  #print(head(intercepts))

  splitter_expression <- paste0('vapply(condition_group, function(z) strsplit(z, "__", fixed=TRUE)[[1]][',
                                seq_along(conditions), '], "a")') %>%
    setNames(conditions)

  #print(splitter_expression)

  x$model_parameters <- x$model_parameters %>%
    dplyr::ungroup() %>%
    dplyr::filter(term != paste0('condition_group', negative_level)) %>%
    dplyr::left_join(intercepts) %>%
    dplyr::left_join(x$mean_variance_relationship) %>%
    dplyr::mutate_(.dots=splitter_expression)

  return ( x )

}

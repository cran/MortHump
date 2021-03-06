% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.codhump.r
\name{summary.codhump}
\alias{summary.codhump}
\alias{print.summary.codhump}
\title{Summary method for codhump models}
\usage{
\method{summary}{codhump}(object, ...)

\method{print}{summary.codhump}(x, ...)
}
\arguments{
\item{object}{the result of a \code{codhump} fit}

\item{...}{further arguments passed to or from other methods}

\item{x}{an object of class "summary.codhump", usually, a result of a call to summary.codhump}
}
\value{
\describe{
 \item{loss}{vector containing the losses in life expectancy due to the hump, first all causes taken together, then for each cause of death.}
 \item{mode}{vector containing the modal age at death from the hump, first all causes taken together, then for each cause of death}
 \item{mean}{vector containing the mean age at death from the hump, first all causes taken together, then for each cause of death}
 \item{sd}{vector containing the standard deviation of the age at death from the hump, first all causes taken together, then for each cause of death}
 }
}
\description{
Produces a summary of a fitted codhump model
}
\details{
This summary function applies to a \code{codhump} object, resulting from a call to the \code{codhump} function.
Using the fitted age- and cause-specific contributions to the hump, it computes summary measures about the intensity, centrality and dispersion of the hump and its components.

\bold{intensity}

The indensity of the hump is measured in terms of the loss of life expectancy implied by the presence of the all-cause hump, as well as each cause-specific contributions.
For the all-cause mortality, it is defined as the difference between life expectancy computed on the fitted force of mortality (sum of both components), and life expectancy computed on the hump-free force of mortality (only the senescence component).
This difference can be interpreted as the mean years of life lost in the population due to the presence of the hump, or alternatively as the potential gain in life expectancy that could be reached in the absence of the hump.
This measure is extended to cause-specific contributions to the hump by substracting separately each of the hump components from the overall force of mortality and recomputing the life expectancy.
This counterfactual approach assumes implicitely the independance of the cause-specific contributions to the hump. In other words, it assumes that removing one cause's hump component will not affect the other causes' components.
In reality, one might argue that people who are under higher risk of death from one cause during early adulthood are also under a higher risk from other causes, and thus that the assumption of independence does not hold.
The loss of life expectancy should thus be interpreted with care and be considered as a purely statistical and theoretical exercice, and not be interpreted as a prediction of what could be gained if we could shelter people from the young adult excess mortality.

\bold{centrality}

Two measures of centrality are currently proposed: the modal and mean age at death from the hump for each cause. By considering the hump as a density, standard summary measures can be applied to it.
The mode is defined as the age at which is observed the highest number of deaths due to the hump, and is therefore not sensitive to the proportion of deaths recorded at other ages.
The mean is defined as the mean age at death for the people who die because of the hump, and is therefore sensitive to the general shape of the hump.
The difference between the mode and the mean can inform on the symetry of the hump, e.g. if the mode is located before the mean this suggests the presence of a flatter hump after the peak.

\bold{dispersion}

Dispersion is measured by the standard deviation of the age at which people die because of the hump.
The all-cause and the cause-specific dispersions are independant of each other and should be interpreted together with measures of centrality.
One may for instance observe a high all-cause dispersion combined with small cause-specific dispersions if the cause-specific humps are not centered on the same age.
Alternatively, if all the cause-specific humps are centered on the same age, the overall dispersion may be smaller than each of the cause-specific dispersions.
Missing values may be caused by remaining negative contributions to the hump, which can be solved by changing parameters in the fit of the \code{\link{codhump}}.
}

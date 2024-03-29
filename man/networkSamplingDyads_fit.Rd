% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R6Class-networkSampling_fit.R
\name{networkSamplingDyads_fit}
\alias{networkSamplingDyads_fit}
\title{Virtual class used to define a family of networkSamplingDyads_fit}
\description{
Virtual class used to define a family of networkSamplingDyads_fit

Virtual class used to define a family of networkSamplingDyads_fit
}
\section{Super class}{
\code{\link[missSBM:networkSampling]{missSBM::networkSampling}} -> \code{networkSamplingDyads_fit}
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{penalty}}{double, value of the penalty term in ICL}

\item{\code{log_lambda}}{double, term for adjusting the imputation step which depends on the type of sampling}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-networkSamplingDyads_fit-new}{\code{networkSamplingDyads_fit$new()}}
\item \href{#method-networkSamplingDyads_fit-show}{\code{networkSamplingDyads_fit$show()}}
\item \href{#method-networkSamplingDyads_fit-update_parameters}{\code{networkSamplingDyads_fit$update_parameters()}}
\item \href{#method-networkSamplingDyads_fit-update_imputation}{\code{networkSamplingDyads_fit$update_imputation()}}
\item \href{#method-networkSamplingDyads_fit-clone}{\code{networkSamplingDyads_fit$clone()}}
}
}
\if{html}{\out{
<details open><summary>Inherited methods</summary>
<ul>
<li><span class="pkg-link" data-pkg="missSBM" data-topic="networkSampling" data-id="print"><a href='../../missSBM/html/networkSampling.html#method-networkSampling-print'><code>missSBM::networkSampling$print()</code></a></span></li>
</ul>
</details>
}}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-networkSamplingDyads_fit-new"></a>}}
\if{latex}{\out{\hypertarget{method-networkSamplingDyads_fit-new}{}}}
\subsection{Method \code{new()}}{
constructor for networkSampling_fit
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{networkSamplingDyads_fit$new(partlyObservedNetwork, name)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{partlyObservedNetwork}}{a object with class partlyObservedNetwork representing the observed data with possibly missing entries}

\item{\code{name}}{a character for the name of sampling to fit on the partlyObservedNetwork}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-networkSamplingDyads_fit-show"></a>}}
\if{latex}{\out{\hypertarget{method-networkSamplingDyads_fit-show}{}}}
\subsection{Method \code{show()}}{
show method
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{networkSamplingDyads_fit$show()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-networkSamplingDyads_fit-update_parameters"></a>}}
\if{latex}{\out{\hypertarget{method-networkSamplingDyads_fit-update_parameters}{}}}
\subsection{Method \code{update_parameters()}}{
a method to update the estimation of the parameters. By default, nothing to do (corresponds to MAR sampling)
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{networkSamplingDyads_fit$update_parameters(...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{use for compatibility}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-networkSamplingDyads_fit-update_imputation"></a>}}
\if{latex}{\out{\hypertarget{method-networkSamplingDyads_fit-update_imputation}{}}}
\subsection{Method \code{update_imputation()}}{
a method to update the imputation of the missing entries.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{networkSamplingDyads_fit$update_imputation(nu)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{nu}}{the matrix of (uncorrected) imputation for missing entries}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-networkSamplingDyads_fit-clone"></a>}}
\if{latex}{\out{\hypertarget{method-networkSamplingDyads_fit-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{networkSamplingDyads_fit$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}

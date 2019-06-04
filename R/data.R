#' American College football
#'
#' Network of American football games between Division IA colleges during regular season Fall 2000. 
#' @docType data
#'
#' @usage data(Football)
#'
#' @format An object of class \code{"list"} with two fields; 
#'  \describe{
#'   \item{X}{network adjacency matrix as a \cite{sparseMatrix} of size 115x115}
#'   \item{label}{vector of teams conferences of size 115 with the following encoding (0 = Atlantic Coast,
#'   1 = Big East,
#'   2 = Big Ten,
#'   3 = Big Twelve,
#'   4 = Conference USA,
#'   5 = Independents,
#'   6 = Mid-American,
#'   7 = Mountain West,
#'   8 = Pacific Ten,
#'   9 = Southeastern,
#'   10 = Sun Belt,
#'   11 = Western Athletic)}
#' }
#'
#' @keywords datasets
#'
#' @references M. Girvan and M. E. J. Newman, Community structure in social and biological networks, Proc. Natl. Acad. Sci. USA 99, 7821-7826 (2002)
#' (\href{https://www.pnas.org/content/pnas/99/12/7821.full.pdf}{PNAS}).
#'
#' @source \href{http://www-personal.umich.edu/~mejn/netdata/}{M. E. J. Newman Network datasets}
#'
#' @examples
#' data(Football)
"Football"



#' Political blogs 
#' 
#' A directed network of hyperlinks between weblogs on US politics, recorded in 2005 by Adamic and Glance. 
#' Only the biggest connected component of the original graph is provided. 
#' 
#' @docType data
#'
#' @usage data(Blogs)
#'
#' @format An object of class \code{"list"} with two fields; 
#'  \describe{
#'   \item{X}{network adjacency matrix as a \cite{sparseMatrix} of size 1222x1222}
#'   \item{label}{vector of political leaning of each blogs (size 1222) with the following encoding (0 = left or liberal,1 = right or conservative)}
#' }
#'
#' @keywords datasets
#'
#' @references Lada A. Adamic and Natalie Glance, "The political blogosphere and the 2004 US Election", in Proceedings of the WWW-2005 Workshop 
#' on the Weblogging Ecosystem (2005) (\href{https://dl.acm.org/citation.cfm?id=1134277}{ACM}).
#'
#' @source \href{http://www-personal.umich.edu/~mejn/netdata/}{M. E. J. Newman Network datasets}
#'
#' @examples
#' data(Blogs)
"Blogs"



#' Books about US politics
#' 
#' A network of books about US politics published around the time of the 2004 presidential election and sold by the online bookseller Amazon.com. Edges between books represent frequent copurchasing of books by the same buyers. 
#' The network was compiled by V. Krebs and is unpublished, but can found on Krebs' web site. Thanks to Valdis Krebs for permission to post these data on this web site.
#' 
#' @docType data
#'
#' @usage data(Books)
#'
#' @format An object of class \code{"list"} with two fields; 
#'  \describe{
#'   \item{X}{network adjacency matrix as a \cite{sparseMatrix} of size 105x105}
#'   \item{label}{ a factor of length  (size 105) with levels "l", "n", or "c" to indicate whether the books are liberal, neutral, or conservative}
#' }
#'
#' @keywords datasets
#'
#'
#' @source \href{http://www-personal.umich.edu/~mejn/netdata/}{M. E. J. Newman Network datasets}
#'
#' @examples
#' data(Books)
"Books"


#' Jazz musicians network
#' 
#' List of edges of the network of Jazz musicians.
#'  
#' @docType data
#'
#' @usage data(Jazz)
#'
#' @format An object of class \code{"sparseMatrix"} with the network adjacency matrix. 
#'
#' @keywords datasets
#'
#' @references P.Gleiser and L. Danon , Community Structure in jazz, Adv. Complex Syst.6, 565 (2003) (\href{https://arxiv.org/abs/cond-mat/0307434}{Arxiv})
#' 
#' @source \href{http://deim.urv.cat/~alexandre.arenas/data/welcome.htm}{A. Arena Network datasets}
#'
#' @examples
#' data(Jazz)
"Jazz"


#' French Parliament votes
#' 
#' 
#'  
#' 
#' @docType data
#'
#' @usage data(FrenchParliament)
#'
#' @format An object of class \code{"list"} with two fields; 
#'  \describe{
#'   \item{X}{matrix of deputy votes a \cite{sparseMatrix} of size 593x570}
#'   \item{labels}{a data frame of deputy meta-data}
#' }
#' @examples
#' data(FrenchParliament)
"FrenchParliament"



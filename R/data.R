#' Ndrangheta mafia covert network dataset
#'
#' Network of co-attendance occurrence attendance of suspected members of the Ndrangheta criminal organization at summits (meetings whose purpose is to make important decisions and/or affiliations, but also to solve internal problems and to establish roles and powers) taking place between 2007 and 2009.

#' @docType data
#'
#' @usage data(Ndrangheta)
#'
#' @format An object of class \code{list} with two fields;
#'  \describe{
#'   \item{X}{network adjacency matrix as a \code{matrix} of size 146x146}
#'   \item{node_meta}{data frame of nodes meta information with features :}
#'   \describe{
#'   \item{Id}{id of the node, rownames of network adjacency matrix}
#'   \item{Locale}{factor with the locali affiliation of the node , "OUT": Suspects not belonging to La Lombardia, "MISS": Information not available, other Locali Id.}
#'   \item{Role}{factor with the type of hierarchical position of the node "MISS": Information not available,"boss": high hierarchical position, "aff": affiliate}
#'  }
#' }
#'
#' @keywords datasets
#'
#' @references{Extended Stochastic Block Models with Application to Criminal Networks, Sirio Legramanti and Tommaso Rigon and Daniele Durante and David B. Dunson, 2021,
#' (\href{https://arxiv.org/abs/2007.08569v2}{arXiv:2007.08569}).}
#'
#' @source \href{https://sites.google.com/site/ucinetsoftware/datasets/covert-networks}{ucinetsoftware/datasets/covert-networks}
#'
#' @examples
#' data(Ndrangheta)
"Ndrangheta"



#' American College football network dataset
#'
#' Network of American football games between Division IA colleges during regular season Fall 2000.
#' @docType data
#'
#' @usage data(Football)
#'
#' @format An object of class \code{list} with two fields;
#'  \describe{
#'   \item{X}{network adjacency matrix as a \code{\link[Matrix]{sparseMatrix}} of size 115x115}
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
#' (\href{https://www.pnas.org/doi/10.1073/pnas.122653799}{PNAS}).
#'
#' @source \href{http://www-personal.umich.edu/~mejn/netdata/}{M. E. J. Newman Network datasets}
#'
#' @examples
#' data(Football)
"Football"


#' Fashion mnist dataset
#'
#' Zalando fashionmnist dataset, sample of 1 000 Zalando's article images from the test set.
#' @docType data
#'
#' @usage data(fashion)
#'
#' @keywords datasets
#'
#' @format An object of class \code{matrix} with a random sample of 1000 images (one per rows) extracted from the fashionmnist dataset.
#'
#'
#' @references Fashion-MNIST: a Novel Image Dataset for Benchmarking Machine Learning Algorithms. Han Xiao, Kashif Rasul, Roland Vollgraf (2017)
#' (\href{https://arxiv.org/abs/1708.07747}{arXiv:1708.07747}).
#'
#' @source \href{https://github.com/zalandoresearch/fashion-mnist}{https://github.com/zalandoresearch/fashion-mnist}
#'
#' @examples
#' data(fashion)
"fashion"







#' Books about US politics network dataset
#'
#' A network of books about US politics published around the time of the 2004 presidential election and sold by the online bookseller Amazon.com. Edges between books represent frequent co-purchasing of books by the same buyers.
#' The network was compiled by V. Krebs and is unpublished, but can found on Krebs' web site. Thanks to Valdis Krebs for permission to post these data on this web site.
#'
#' @docType data
#'
#' @usage data(Books)
#'
#' @format An object of class \code{list} with two fields;
#'  \describe{
#'   \item{X}{network adjacency matrix as a \code{\link[Matrix]{sparseMatrix}} of size 105x105}
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


#' Jazz musicians network dataset
#'
#' List of edges of the network of Jazz musicians.
#'
#' @docType data
#'
#' @usage data(Jazz)
#'
#' @format An object of class \code{\link[Matrix]{sparseMatrix}} with the network adjacency matrix.
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




#' Mushroom data
#'
#' Categorical data from UCI Machine Learning Repository describing 8124
#' mushrooms with 22 phenotype variables. Each mushroom is classified as "edible"
#' or "poisonous" and the goal is to recover the mushroom class from its
#' phenotype.
#'
#'
#' @docType data
#'
#' @keywords datasets
#'
#' @usage data(mushroom)
#'
#' @format An R data.frame with a variable edibility used as label and 22
#'  categorical variables with no names. More detail on the UCI webpage
#'  describing the data.
#'
#' @source \url{https://archive.ics.uci.edu/ml/datasets/Mushroom}
#' @examples
#' data(mushroom)
"mushroom"


#' Young People survey data
#'
#' Young people survey data from Miroslav Sabo and available on the Kaggle
#' platform. This is an authentic example of questionnaire data where Slovakian
#' young people (15-30 years old) were asked musical preferences according to
#' different genres (rock, hip-hop, classical, etc.).
#'
#' @docType data
#'
#' @usage data(Youngpeoplesurvey)
#'
#' @format An R data.frame with columns containing each of the 150 original
#'   variables of the study.
#'
#'
#' @keywords datasets
#'
#' @source \url{https://www.kaggle.com/miroslavsabo/young-people-survey}
#'
#' @examples
#' data(Youngpeoplesurvey)
"Youngpeoplesurvey"




#' Fifa data
#'
#' A random sample of 6000 players from the FIFA videogame with various statistics on all player ranging
#' from position, cost in the game, capacity in offense/defense, speed, etc. 
#' Two columns pos_x, pos_y with average player possible positions (in opta coordiantes)
#' were derived from the raw data.  was also u.
#'
#' @docType data
#'
#' @usage data(Fifa)
#'
#' @format An R data.frame with columns containing each of the descriptive
#'   statistics of a player.
#'
#'
#' @keywords datasets
#'
#' @source \url{https://www.kaggle.com/stefanoleone992/fifa-20-complete-player-dataset?select=players_20.csv}
#'
#' @examples
#' data(Fifa)
"Fifa"





#' NewGuinea data
#'
#' \code{\link{NewGuinea}} a social network of 16 tribes, where two types of interactions were recorded, amounting to either friendship or enmity [read-cultures-1954]. 
#'
#' @docType data
#'
#' @usage data(NewGuinea)
#'
#' @format A binary array of size (16,16,3) the first layer encodes enmity, the second, the friendship relations. The third, no relations between the two tribes.
#' @references Kenneth E. Read, “Cultures of the Central Highlands, New Guinea”, Southwestern J. of Anthropology, 10(1):1-43 (1954). DOI: 10.1086/soutjanth.10.1.3629074
#'
#' @keywords datasets
#'
#' @source \url{https://networks.skewed.de/net/new_guinea_tribes}
#'
#' @examples
#' data(NewGuinea)
"NewGuinea"

#' SevenGraders data
#'
#' \code{\link{SevenGraders}} A small multiplex network of friendships among 29 seventh grade students in Victoria, Australia. Students nominated classmates for three different activities (who do you get on with in the class, who are your best friends, and who would you prefer to work with). Edge direction for each of these three types of edges indicates if node i nominated node j, and the edge weight gives the frequency of this nomination. Students 1-12 are boys and 13-29 are girls. The KONECT version of this network is the collapse of de Domenico's multiplex version.
#' @docType data
#'
#' @usage data(SevenGraders)
#'
#' @format A binary array of size (29,29,3) containing directed graphs. The first layer encodes "getting along in class" while the second encodes the best-friendship (can be one-way). The third encodes the preferred work relation.
#'
#' @keywords datasets
#'
#' @source \url{https://networks.skewed.de/net/7th_graders}
#' @references M. Vickers and S. Chan, "Representing Classroom Social Structure." Melbourne: Victoria Institute of Secondary Education, (1981).
#' @examples
#' data(SevenGraders)
"SevenGraders"



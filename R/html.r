#' Get HTML href content
#'
#' @param file character. Path to HTML file.
#' @param css character. HTML elements (one of css) to select. Default \code{'body a'}
#'
#' @importFrom rvest read_html html_elements html_attr
#'
#' @return A character vector.
#' @export
#'
#' @examples
#' \donttest{
#' Download('https://ftp.ensembl.org/pub', 'test.html')
#' content = htmlHref('test.html')
#' head(content)
#' }
#' 
htmlHref = function(file, css = 'body a') 
  rvest::html_attr(rvest::html_elements(rvest::read_html(file), css), 'href')

#' Get HTML text content
#'
#' @param file character. Path to HTML file.
#' @param css character. HTML elements (one of css) to select. Default \code{'body'}
#'
#' @importFrom rvest read_html html_elements html_text
#' 
#' @return A character vector.
#' @export
#'
#' @examples
#' \donttest{
#' Download('http://www.virusite.org', 'test.html')
#' content = htmlText('test.html')
#' head(content)
#' }
#' 
htmlText = function(file, css = 'body') 
  rvest::html_text(rvest::html_elements(rvest::read_html(file), css))


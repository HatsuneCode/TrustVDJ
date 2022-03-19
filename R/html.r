#' Get HTML href content
#'
#' @param file character. Path to HTML file.
#' @param css character. HTML elements (one of css) to select. Default \code{'body a'}
#'
#' @importFrom rvest read_html html_elements html_text
#'
#' @return
#' @export
#'
#' @examples
#' \donttest{
#' Download('https://ftp.ensembl.org/pub', 'test.html')
#' content = htmlText('test.html', css = 'body a')
#' head(content)
#' }
#' 
htmlHref = function(file, css = 'body a') rvest::html_attr(rvest::html_elements(rvest::read_html(file), css), 'href')

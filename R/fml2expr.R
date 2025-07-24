#' Turn sum fromula vector into expression vector
#'
#' @param x Character vector containing chemical formulas.
#' @param expr Will return vector of expressions when TRUE.
#'
#' @returns
#' @export
#'
#' @examples
#' text <- c("Beispiele: A1, Bc23, D456, Efg7890, Z9A12.", "C6H12O6","[M+H]+")
#' fml2expr(text)
#' fml2expr(gsub(" ", "~", c("Text C6H12Cl6 und H3Na text", "A1B2", "Kein Match hier")), expr=T)
#' texts <- c(
#'   "C6H12Cl6",         # einfache Formel
#'   "H3Na",             # einfache Formel
#'   "Text mit Leerzeichen",
#'   "A1B2",             # direkt aufeinanderfolgende Elemente
#'   "Fe(NO3)3",         # Gruppe mit Index
#'   "SO4^2-",           # Ladung mit Oxidationszahl
#'   "Na+", "Cl-",       # einfache Ionen
#'   "H₂O",              # Unicode-Subscript
#'   "Ca(OH)2"           # Gruppe mit Index
#' )
#' 
#' # Als Ausdruck für Plot
#' sapply(texts , function(x) { inherits(try(fml2expr(x, expr = TRUE), silent=TRUE), "try-error") })
#' exprs <- fml2expr(texts, expr = TRUE)
#' 
#' # Plot mit Legende
#' plot(1:10, 1:10, pch = 16, col = 1:10, xlim = c(0, 11), ylim = c(0, 11))
#' legend("topright", legend = exprs, col = 1:10, pch = 16, title = "Formeln")
#' text(2,2,labels=fml2expr(texts[2]))
#' text(3:4,3:4,labels=fml2expr(texts[3:4]))

fml2expr <- function(x, expr = TRUE) {
  unicode_subs <- c("₀"="0", "₁"="1", "₂"="2", "₃"="3", "₄"="4",
                    "₅"="5", "₆"="6", "₇"="7", "₈"="8", "₉"="9")
  
  process_element <- function(text) {
    # 1. Unicode-Subscripts ersetzen
    for (u in names(unicode_subs)) {
      text <- gsub(u, unicode_subs[[u]], text, fixed = TRUE)
    }
    
    # 2. Leerzeichen durch ~ ersetzen
    text <- gsub(" ", "~", text)
    
    # 3. Gruppen in Klammern mit Index: (OH)2 → (OH)[2]
    text <- gsub("\\(([A-Za-z0-9]+)\\)(\\d{1,3})", "(\\1)[\\2]", text, perl = TRUE)
    
    # 4. Chemische Formeln: Großbuchstabe + optional Kleinbuchstabe + Zahl → A[1]
    text <- gsub("([A-Z][a-z]?)(\\d{1,3})", "\\1[\\2]", text, perl = TRUE)
    
    # 5. Ladungen:
    # a) Oxidationszahlen wie ^2-, ^3+ → ^'2-', ^'3+'
    text <- gsub("\\^(\\d+[+-])", "^'\\1'", text, perl = TRUE)
    
    # b) einfache + oder - am Ende → ^'+', ^'-'
    text <- gsub("([A-Za-z0-9\\]\\)])([+-])$", "\\1^'\\2'", text, perl = TRUE)
    
    # c) einfache + oder - innerhalb → ^'+', ^'-'
    #text <- gsub("([A-Za-z0-9\\]\\)])([+-])(?![0-9])", "\\1^'\\2'", text, perl = TRUE)
    
    # 6. Füge ~ zwischen benachbarten chemischen Einheiten ein:
    text <- gsub("\\](?=\\()", "]~", text, perl = TRUE)
    text <- gsub("\\](?=[A-Z])", "]~", text, perl = TRUE)
    
    # 7. Entferne überflüssiges ~ am Ende
    text <- sub("~$", "", text)
    
    # Rückgabe
    if (expr) {
      parse(text = text)[[1]]
    } else {
      text
    }
  }
  n <- length(x)
  if (expr) {
    if (n==1) process_element(x) else lapply(x, process_element)
  } else {
    if (n==1) process_element(x) else vapply(x, process_element, FUN.VALUE = character(1))
  }
}

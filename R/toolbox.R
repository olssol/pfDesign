#' Export results into template latex tables
#'
#'
#' @export
#'
tkExpRst <- function(numbers, template_f, out_f, sub.str="AA", append = FALSE) {
    if (!file.exists(template_f)) {
        return;
    }
    ##read
    tpla <- readChar(template_f, file.info(template_f)$size);

    ##substitute
    for (i in seq_len(length(numbers))) {
        tpla <- sub(sub.str, numbers[i], tpla);
    }

    ##write out
    write(tpla, file = out_f, append = append)
}

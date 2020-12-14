.check_file <- function(filename, fileType = "File") {
  if (!file.exists(filename))
    stop(fileType, " ", filename, " does not exist.", call. = FALSE)
  filename <- normalizePath(filename)
  return(filename)
}

.check_integer <- function(int, arg_name) {
  if (!(is.numeric(int)))
    stop("`", arg_name, "` must be a positive integer.", call. = FALSE)
  if (any(int <= 0) | !all(as.integer(int) == int))
    stop("`", arg_name, "` must be a positive integer.", call. = FALSE)
  return(int)
}  

.check_double <- function(d, arg_name) {
  if (!is.numeric(d))
    stop("`", arg_name, "` must be a double.", call. = FALSE)
  return(d)
}

.check_pos_double <- function(d, arg_name) {
  if (!is.numeric(d) | any(d <= 0))
    stop("`", arg_name, "` must be a positive double.", call. = FALSE)
  return(d)
}

.check_range <- function(r, arg_name) {
  if (length(r) != 2)
    stop("`", arg_name, "` must be a vector of size 2.", call. = FALSE)
  return(r)
}

.check_char <- function(str, arg_name) {
  if (!is.character(str)) {
    stop("`", arg_name, "` must be a string.", call. = FALSE)
  }
  return(str)
}

.check_bool <- function(bool, arg_name) {
  if (!is.logical(bool))
    stop("`", arg_name, "` must be logical.", call. = FALSE)
  return(bool)
}

.check_in_vector <- function(arg, arg_name, allowed_args, type = "error") {
  if (!all(is.element(arg, allowed_args))) {
    if (type == "error") {
      stop("Invalid `", arg_name, "`. Available elements: \n  - ",
           paste(allowed_args, collapse="\n  - "),
           call. = FALSE)
    } else if (type == "warning") {
      warning("Invalid `", arg_name, "`. Available elements: \n  - ",
              paste(allowed_args, collapse="\n  - "),
              call. = FALSE)
    }
  }
  arg <- arg[arg %in% allowed_args]
  return(arg)
}

.check_object_class <- function(arg, arg_name, arg_class) {
  if (any(class(arg) != arg_class))
    stop("`", arg_name, "` must be a `", arg_class, "` object.")
  return(arg)
}
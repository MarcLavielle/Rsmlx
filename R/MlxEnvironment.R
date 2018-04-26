MlxEnvironment <- new.env()

assign("MLX_DIRECTORY", "", envir = MlxEnvironment)
assign("MLXCONNECTORS_LIB_NAME", "", envir = MlxEnvironment)
assign("MLXCONNECTORS_LIB_PATH", "", envir = MlxEnvironment)

OS = .getOS()
if (OS == "Unix") {
  assign("SYSTEM_PATH", Sys.getenv("LD_LIBRARY_PATH"), envir = MlxEnvironment)
} else if (OS == "Apple") {
  # path is set during installation
} else if (OS == "Windows") {
  assign("SYSTEM_PATH", Sys.getenv("PATH"), envir = MlxEnvironment)
}
remove(OS)
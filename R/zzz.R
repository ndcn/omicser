.onLoad <- function(libname = find.package("omicser"), pkgname = "omicser") {
  # CRAN Notes avoidance
  if(getRversion() >= "2.15.1") {
    utils::globalVariables(
      c(
        # gen_config_table
        "UI",
        "field",
        # mod_pg_diff_expr_server
        "test_type",
        "versus",
        "f",
        "significant",
        "point_color",
        "neglogpval",
        "ID",
        "UI",
        # mod_pg_expression_server
        "grp"
      )
    )
  }
  invisible()
}

# ELM API tests — all HTTP calls use vcr cassettes so no real network needed.
# To record fresh cassettes: set record = "new_episodes" in helper-vcr.R,
# then run devtools::test() with internet access.

test_that("query_elm returns NULL gracefully when cassette returns error", {
  # Use webmockr to stub a failing request without a cassette file
  webmockr::enable()
  on.exit(webmockr::disable(), add = TRUE)

  webmockr::stub_request("post", uri = "https://elm.eu.org/search/sequence/") %>%
    webmockr::to_return(status = 503, body = "Service Unavailable")

  result <- query_elm("MDFDY", taxon = "Viruses", timeout = 5L)
  expect_null(result)
})

test_that("query_elm_cached returns same result on repeated calls (no extra HTTP)", {
  webmockr::enable()
  on.exit({
    webmockr::disable()
    # Clear the memoise cache between tests
    memoise::forget(query_elm_cached)
  }, add = TRUE)

  fake_json <- '{"instances": [], "classes": []}'
  webmockr::stub_request("post", uri = "https://elm.eu.org/search/sequence/") %>%
    webmockr::to_return(status = 200, body = fake_json,
                        headers = list("Content-Type" = "application/json"))

  r1 <- query_elm_cached("MDFDY", taxon = "Viruses", timeout = 5L)
  r2 <- query_elm_cached("MDFDY", taxon = "Viruses", timeout = 5L)

  # Both calls should return same object; second call hits cache, not stub
  expect_identical(r1, r2)
  # Stub was only hit once (memoise served the second)
  expect_equal(length(webmockr::request_registry()$request_signatures$stack), 1L)
})

test_that("build_concat produces correct positions", {
  frags <- c(A = "MDFDY", B = "LKPTG")
  result <- build_concat(frags)
  expect_named(result, c("concat", "map"))
  expect_equal(result$map$isolate, c("A", "B"))
  # Start of A is 1; end = nchar("MDFDY") = 5
  expect_equal(result$map$start[1], 1L)
  expect_equal(result$map$end[1], 5L)
  # B starts after linker
  linker_len <- nchar(potytools:::LINKER)
  expect_equal(result$map$start[2], 5L + linker_len + 1L)
})

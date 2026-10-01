// Regression tests for reproducible quantitative-validation CLI execution.

test "cli_solver_filter_controls_actual_execution" {
  assert_true(cli_solver_selected("all", "chorin"))
  assert_true(cli_solver_selected("all", "simple"))
  assert_true(cli_solver_selected("chorin", "chorin"))
  assert_false(cli_solver_selected("chorin", "simple"))
  assert_false(cli_solver_selected("simple", "pcg"))
}

test "checkpoint_steps_must_be_positive_and_strictly_increasing" {
  let checkpoints = parse_checkpoint_steps("2,5,10")
  assert_eq(checkpoints.length(), 3)
  assert_eq(checkpoints[0], 2)
  assert_eq(checkpoints[1], 5)
  assert_eq(checkpoints[2], 10)
  assert_eq(parse_checkpoint_steps("10,5").length(), 0)
  assert_eq(parse_checkpoint_steps("5,5").length(), 0)
  assert_eq(parse_checkpoint_steps("0,5").length(), 0)
  assert_eq(parse_checkpoint_steps("bad").length(), 0)
}

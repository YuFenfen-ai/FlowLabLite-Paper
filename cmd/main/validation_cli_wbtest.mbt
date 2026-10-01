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

test "chorin_runtime_dt_can_be_changed_for_grid_validation" {
  let original = get_dt()
  set_chorin_dt(0.0005)
  assert_eq(get_dt(), 0.0005)
  set_chorin_dt(original)
}

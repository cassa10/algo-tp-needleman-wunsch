def print_test(n):
    print(f"[Test {n}]")


def print_test_pass(n, is_passed):
    if is_passed:
        print(f"[Test {n}] - PASSED")
    else:
        print(f"[Test {n}] - ERROR")
    print(f"------------------------------")


def execute_tests(tests):
    for i, test in enumerate(tests):
        n = i + 1
        print_test(n)
        print_test_pass(n, test())

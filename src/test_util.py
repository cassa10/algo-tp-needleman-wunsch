def print_test(n):
    print(f"[Test {n}]")


def print_test_pass(n, is_passed):
    if is_passed:
        print(f"[Test {n}] - PASSED")
    else:
        print(f"[Test {n}] - ERROR")
    print(f"------------------------------")


def print_test_all_passed(is_all_passed):
    if is_all_passed:
        print(f"All tests passed")
    else:
        print("ERROR - No all tests passed")


def execute_tests(tests):
    all_passed = True
    for i, test in enumerate(tests):
        n = i + 1
        print_test(n)
        passed = test()
        all_passed = all_passed and passed
        print_test_pass(n, passed)
    print_test_all_passed(all_passed)

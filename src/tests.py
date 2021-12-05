from src import needleman_wunsh_test as nw_test, msa_test


def run_all():
    all_passed = True
    all_passed = all_passed and nw_test.run_tests()
    return all_passed and msa_test.run_tests()

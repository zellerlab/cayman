from . import test_annotation_cli


def load_tests(loader, suite, pattern):
    suite.addTests(loader.loadTestsFromModule(test_annotation_cli))
    return suite
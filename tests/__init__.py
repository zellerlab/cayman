from . import test_annotation_cli, test_hmmloader


def load_tests(loader, suite, pattern):
    suite.addTests(loader.loadTestsFromModule(test_hmmloader))
    suite.addTests(loader.loadTestsFromModule(test_annotation_cli))
    return suite
import os
import unittest

from tests.test_config import BaseConfig
from lib.api import omim
from lib.model import config


class OmimTest(BaseConfig):

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.omim = omim.Omim(mimdir=cls.input_path)

    def test_init(self):
        omim_obj = omim.Omim(
            mimdir=self.input_path
        )
        self.assertIsNotNone(omim_obj)

    def test_gene_omim_to_entrez(self):
        tests = [
            ("603024", "8289"),  # ARID1A
            ("600014", "6595"),  # SMARCA2
        ]
        for i, (test, correct) in enumerate(tests):
            with self.subTest(i=i):
                self.assertEqual(
                    self.omim.mim_gene_to_entrez_id(test),
                    correct
                )

    def test_entrez_to_symbol(self):
        tests = [
            ("8289", "ARID1A"),  # ARID1A
            ("6595", "SMARCA2"),  # SMARCA2
        ]
        for i, (test, correct) in enumerate(tests):
            with self.subTest(i=i):
                self.assertEqual(
                    self.omim.entrez_id_to_symbol(test),
                    correct
                )

    def test_pheno_to_genes(self):
        tests = [
            ("135900", {
                "614556": {
                    "gene_id": "57492",
                    "gene_symbol": "ARID1B",
                    "gene_omim_id": "614556"
                }
            }),
        ]
        for i, (test, correct) in enumerate(tests):
            with self.subTest(i=i):
                self.assertDictEqual(
                    self.omim.mim_pheno_to_gene(test),
                    correct
                )

    def test_pheno_to_ps(self):
        tests = [
            ("135900", "PS135900"),
            ("614607", "PS135900")
        ]
        for i, (test, correct) in enumerate(tests):
            with self.subTest(i=i):
                result = self.omim.omim_id_to_phenotypic_series(test)
                self.assertEqual(result, correct)

    def test_deprecated_ids(self):
        tests = [
            ("106160", ["106165"]),
            (["106200"], ["106210"]),
            ("107253", []),
            (["106200", "117200"], ["106210"]),
            ("135900", ["135900"]),
        ]

        for test, correct in tests:
            with self.subTest(i=test):
                replaced = self.omim.replace_deprecated_all(test)
                self.assertListEqual(replaced, correct)

    def test_pheno_to_name(self):
        tests = [
            ("135900", "COFFIN-SIRIS SYNDROME 1; CSS1"),
            ("106200", "MOVED TO 106210"),
        ]
        for test, correct in tests:
            with self.subTest(i=test):
                res = self.omim.mim_pheno_to_syndrome_name(test)
                self.assertEqual(res, correct)

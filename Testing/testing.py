import unittest
import pandas as pd
import functions


class MyTestCase(unittest.TestCase):
    def __init__(self, test_data):
        super().__init__()
        self.test_data = test_data

    def tranform_data_test(self):
        output_data = tranform_data(self.test_data) # hier richtige funktion einfügen

        expected_data = pd.DataFrame({})

        pd.testing.assert_frame_equal(output_data, expected_data)
        pass
        # TODO add transformation test


if __name__ == '__main__':
    unittest.main()

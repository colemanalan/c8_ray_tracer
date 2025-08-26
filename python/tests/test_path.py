import unittest
from c8_tracer import Vec3
from c8_tracer.path import Path, SignalPath


class TestPath(unittest.TestCase):
    def test_single_point_constructor(self):
        p0 = Vec3(1.0, 2.0, 3.0)
        path = Path(p0)
        self.assertEqual(path.getStart(), p0)
        self.assertEqual(path.getEnd(), p0)
        self.assertEqual(len(path), 0)  # No segments yet
        self.assertEqual(path.getLength(), 0.0)

    def test_multiple_point_constructor(self):
        points = [Vec3(0, 0, 0), Vec3(1, 0, 0), Vec3(1, 1, 0)]
        path = Path(points)
        self.assertEqual(path.getStart(), points[0])
        self.assertEqual(path.getEnd(), points[-1])
        self.assertEqual(len(path), 2)
        self.assertEqual(path.getPoint(1), points[1])
        self.assertEqual(list(path), points)

    def test_add_to_end(self):
        path = Path(Vec3(0, 0, 0))
        path.addToEnd(Vec3(1, 0, 0))
        path.addToEnd(Vec3(1, 1, 0))
        self.assertEqual(len(path), 2)
        self.assertEqual(path.getEnd(), Vec3(1, 1, 0))
        self.assertEqual(path.getPoint(2), Vec3(1, 1, 0))

    def test_remove_from_end(self):
        points = [Vec3(0, 0, 0), Vec3(1, 0, 0), Vec3(1, 1, 0)]
        path = Path(points)
        path.removeFromEnd()
        self.assertEqual(len(path), 1)
        self.assertEqual(path.getEnd(), Vec3(1, 0, 0))
        path.removeFromEnd()
        self.assertEqual(len(path), 0)
        self.assertEqual(path.getEnd(), Vec3(0, 0, 0))

    def test_iteration_and_len(self):
        points = [Vec3(0, 0, 0), Vec3(1, 0, 0), Vec3(1, 1, 0)]
        path = Path(points)
        self.assertEqual(len(path), 2)
        self.assertEqual(list(path), points)

    def test_repr(self):
        path = Path(Vec3(0, 0, 0))
        repr_str = repr(path)
        self.assertTrue("Path" in repr_str or "<" in repr_str)  # Basic sanity check


class TestSignalPath(unittest.TestCase):
    def setUp(self):
        self.points = [Vec3(0, 0, 0), Vec3(1, 0, 0), Vec3(1, 1, 0)]
        self.signal_path = SignalPath(
            1.23,
            1.5,
            1.0,
            1.33,
            Vec3(0.1, 0.2, 0.3),
            Vec3(2.0, 2.0, 2.0),
            550.0,
            self.points,
            0.0,
            1.23,
        )

    def test_attributes(self):
        self.assertAlmostEqual(self.signal_path.propagation_time, 1.23)
        self.assertAlmostEqual(self.signal_path.average_refractive_index, 1.5)
        self.assertAlmostEqual(self.signal_path.refractive_index_source, 1.0)
        self.assertAlmostEqual(self.signal_path.refractive_index_destination, 1.33)
        self.assertEqual(self.signal_path.emit, Vec3(0.1, 0.2, 0.3))

    def test_inherited_path_behavior(self):
        self.assertEqual(self.signal_path.getStart(), self.points[0])
        self.assertEqual(self.signal_path.getEnd(), self.points[-1])
        self.assertEqual(len(self.signal_path), 2)
        self.assertEqual(list(self.signal_path), self.points)

    def test_repr(self):
        repr_str = repr(self.signal_path)
        self.assertTrue("SignalPath" in repr_str or "<" in repr_str)

    def test_iteration(self):
        iterated_points = [p for p in self.signal_path]
        self.assertEqual(iterated_points, self.points)

    def test_mutability(self):
        self.signal_path.propagation_time = 2.0
        self.signal_path.average_refractive_index = 1.6
        self.assertEqual(self.signal_path.propagation_time, 2.0)
        self.assertEqual(self.signal_path.average_refractive_index, 1.6)


if __name__ == "__main__":
    unittest.main()

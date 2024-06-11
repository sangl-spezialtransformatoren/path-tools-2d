from typing import Self, Union, Iterable, TypeAlias


class Point:
    x: float
    y: float


type PointLike = Point | (float, float)


class Segment:
    length: float
    start: Point
    end: Point

    def reversed(self) -> Self: ...

    def point_at_t(self, t: float) -> Point: ...

    def t_of_point(self, point: PointLike) -> float: ...

    def slice(self, t1: float, t2: float) -> Self: ...

    def contains(self, point: PointLike) -> bool: ...

    def intersect(self, other: Segment) -> list[Point]: ...


class Arc(Segment):
    rx: float
    ry: float
    phi: float
    large_arc: bool
    sweep: bool

    def __init__(
            self,
            start: PointLike,
            end: PointLike,
            rx: float,
            ry: float,
            phi: float,
            large_arc: bool,
            sweep: bool): ...

    @classmethod
    def from_center_params(
            cls,
            center: PointLike,
            rx: float,
            ry: float,
            phi: float,
            theta1: float,
            d_theta: float
    ) -> Self: ...

    def get_center(self) -> Point: ...


class Line(Segment):
    def __init__(self, start: PointLike, end: PointLike): ...


type AnySegment = Union[Line, Arc]


def combine_segments(segments: Iterable[AnySegment]) -> list[list[AnySegment]]: ...


PointTuple: TypeAlias = tuple[float, float]


def points_to_beziers(points: Iterable[PointTuple], tol: float) \
        -> list[tuple[PointTuple, PointTuple, PointTuple, PointTuple]]: ...


def ramer_douglas_peucker(segments: list[AnySegment], tol: float) -> list[AnySegment]: ...

// Original Work Copyright (c) 2010 Amit J Patel, Modified Work Copyright (c) 2016 Jay M Stevens

#include "PolygonalMapGeneratorPrivatePCH.h"
#include "Maps/PointGenerators/PointGenerator.h"
#include "Maps/PolygonMap.h"
#include "MapDebugVisualizer.h"
#include "Maps/MapDataHelper.h"
#include "Classes/Diagrams/Voronoi.h"

DEFINE_LOG_CATEGORY(LogWorldGen);

UPolygonMap::UPolygonMap()
{
	DefaultCenter = new FMapCenter();
	DefaultCorner = new FMapCorner();
	DefaultEdge = new FMapEdge();
}

UPolygonMap::~UPolygonMap()
{
	delete DefaultCenter;
	delete DefaultCorner;
	delete DefaultEdge;
}

void UPolygonMap::CreatePoints(UPointGenerator* pointSelector, const int32& numberOfPoints)
{
	if (pointSelector == NULL)
	{
		return;
	}

	PointSelector = pointSelector;
	Points = PointSelector->GeneratePoints(numberOfPoints);
	MinPointLocation = PointSelector->MinPoint();
	MaxPointLocation = PointSelector->MaxPoint();
}

void UPolygonMap::BuildGraph(const int32& mapSize, const FWorldSpaceMapData& data)
{
	MapSize = mapSize;
	MapData = data;
	Voronoi voronoi(Points);
	for (int i = 0; i < voronoi.sites.Num(); i++)
	{
		VSite site = voronoi.sites[i];

		// Create centers, corners, and edges
		for (auto siteEdge : site.edges)
		{
			// Voronoi Edge (corner to corner)
			FVector4 vEdge = siteEdge.vEdge;
			// Delaunay Edge (center to center)
			FVector4 dEdge = siteEdge.dEdge;

			// Corners come from edge vertices
			FMapCorner* cornerOne = MakeCorner(FVector2D(vEdge.X, vEdge.Y));
			FMapCorner* cornerTwo = MakeCorner(FVector2D(vEdge.Z, vEdge.W));
			FMapCenter* centerOne = MakeCenter(FVector2D(dEdge.X, dEdge.Y));
			FMapCenter* centerTwo = MakeCenter(FVector2D(dEdge.Z, dEdge.W));

			FMapEdge* edge = new FMapEdge();
			edge->DelaunayEdge0 = centerOne;
			edge->DelaunayEdge1 = centerTwo;
			edge->VoronoiEdge0 = cornerOne;
			edge->VoronoiEdge1 = cornerTwo;
			edge->Midpoint = FVector2D(FMath::Lerp(vEdge.X, vEdge.Z, 0.5f), FMath::Lerp(vEdge.Y, vEdge.W, 0.5f));
			edge->Index = Edges.Num();
			Edges.Add(edge);

			// Corners point to edges
			cornerOne->Protrudes.AddUnique(edge);
			cornerTwo->Protrudes.AddUnique(edge);
			// Centers point to edges
			centerOne->Borders.AddUnique(edge);
			centerTwo->Borders.AddUnique( edge);

			// Centers point to centers
			centerOne->Neighbors.AddUnique(centerTwo);
			centerTwo->Neighbors.AddUnique(centerOne);

			// Corners point to corners
			cornerOne->Adjacent.AddUnique(cornerTwo);
			cornerTwo->Adjacent.AddUnique(cornerOne);

			// Centers point to corners
			centerOne->Corners.AddUnique(cornerOne);
			centerOne->Corners.AddUnique(cornerTwo);
			centerTwo->Corners.AddUnique(cornerOne);
			centerTwo->Corners.AddUnique(cornerTwo);

			// Corners point to centers
			cornerOne->Touches.AddUnique(centerOne);
			cornerOne->Touches.AddUnique(centerTwo);
			cornerTwo->Touches.AddUnique(centerOne);
			cornerTwo->Touches.AddUnique(centerTwo);
		}
	}
	UE_LOG(LogWorldGen, Log, TEXT("Created a total of %d Centers, %d Corners, and %d Edges."), Centers.Num(), Corners.Num(), Edges.Num());
}

FMapCenter* UPolygonMap::MakeCenter(const FVector2D& point)
{
	if (CenterLookup.Contains(point))
	{
		return CenterLookup[point];
	}

	FMapCenter* center = new FMapCenter();
	center->CenterData.Point = point;
	center->CenterData = UMapDataHelper::RemoveBorder(center->CenterData);
	center->CenterData = UMapDataHelper::RemoveOcean(center->CenterData);
	center->Index = Centers.Num();

	Centers.Add(center);
	CenterLookup.Add(point, center);

	return center;
}

FMapCenter& UPolygonMap::GetCenter(const int32& index) const
{
	if (index < 0)
	{
		return *DefaultCenter;
	}
	return *Centers[index];
}

FMapCorner* UPolygonMap::MakeCorner(const FVector2D& point)
{
	if (CornerLookup.Contains(point))
	{
		return CornerLookup[point];
	}

	FMapCorner* corner = new FMapCorner();
	corner->CornerData.Point = point;
	corner->Index = Corners.Num();
	corner->CornerData = UMapDataHelper::RemoveOcean(corner->CornerData);
	if (PointSelector->PointIsOnBorder(corner->CornerData.Point))
	{
		corner->CornerData = UMapDataHelper::SetBorder(corner->CornerData);
	}
	else
	{
		corner->CornerData = UMapDataHelper::RemoveBorder(corner->CornerData);
	}
	
	Corners.Add(corner);
	CornerLookup.Add(point, corner);

	return corner;
}

FMapCorner& UPolygonMap::GetCorner(const int32& index) const
{
	if (index < 0)
	{
		return *DefaultCorner;
	}
	return *Corners[index];
}

FMapEdge& UPolygonMap::GetEdge(const int32& index) const
{
	if (index < 0)
	{
		return *DefaultEdge;
	}
	return *Edges[index];
}
FMapEdge& UPolygonMap::FindEdgeFromCenters(const FMapCenter& v0, const FMapCenter& v1) const
{
	if (v0.Index < 0 || v0.Index >= Edges.Num() || v1.Index <0 || v1.Index >= Edges.Num())
	{
		return GetEdge(-1);
	}
	for (int i = 0; i < Edges.Num(); i++)
	{
		if ((*Edges[i]->DelaunayEdge0 == v0 && *Edges[i]->DelaunayEdge1 == v1) || (*Edges[i]->DelaunayEdge0 == v1 && *Edges[i]->DelaunayEdge1 == v0))
		{
			return GetEdge(i);
		}
	}
	return GetEdge(-1);
}
FMapEdge& UPolygonMap::FindEdgeFromCorners(const FMapCorner& v0, const FMapCorner& v1) const
{
	if (v0.Index < 0 || v0.Index >= Edges.Num() || v1.Index <0 || v1.Index >= Edges.Num())
	{
		return GetEdge(-1);
	}
	for (int i = 0; i < Edges.Num(); i++)
	{
		if ((*Edges[i]->VoronoiEdge0 == v0 && *Edges[i]->VoronoiEdge1 == v1) || (*Edges[i]->VoronoiEdge0 == v1 && *Edges[i]->VoronoiEdge1 == v0))
		{
			return GetEdge(i);
		}
	}
	return GetEdge(-1);
}

FMapCenter& UPolygonMap::FindCenterFromCorners(FMapCorner CornerA, FMapCorner CornerB) const
{
	for (int i = 0; i < CornerA.Touches.Num(); i++)
	{
		for (int j = 0; j < CornerB.Touches.Num(); j++)
		{
			if (CornerA.Touches[i] == CornerB.Touches[j])
			{
				return *CornerA.Touches[i];
			}
		}
	}
	return GetCenter(-1);
}

int32 UPolygonMap::GetCenterNum() const
{
	return Centers.Num();
}
int32 UPolygonMap::GetCornerNum() const
{
	return Corners.Num();
}
int32 UPolygonMap::GetEdgeNum() const
{
	return Edges.Num();
}

int32 UPolygonMap::GetGraphSize() const
{
	return MapSize;
}

/*void UPolygonMap::UpdateCenter(const FMapCenter& center)
{
	if (center.Index < 0 || center.Index >= Centers.Num())
	{
		return;
	}
	Centers[center.Index] = center;
}
void UPolygonMap::UpdateCorner(const FMapCorner& corner)
{
	if (corner.Index < 0 || corner.Index >= Corners.Num())
	{
		return;
	}
	Corners[corner.Index] = corner;
}
void UPolygonMap::UpdateEdge(const FMapEdge& edge)
{
	if (edge.Index < 0 || edge.Index >= Edges.Num())
	{
		return;
	}
	Edges[edge.Index] = edge;
}*/

void UPolygonMap::ImproveCorners()
{
	TArray<FVector2D> newCorners;
	newCorners.SetNumZeroed(Corners.Num());

	// First we compute the average of the centers next to each corner.
	for (int i = 0; i < Corners.Num(); i++)
	{
		if (UMapDataHelper::IsBorder(Corners[i]->CornerData))
		{
			newCorners[i] = Corners[i]->CornerData.Point;
		}
		else
		{
			FVector2D point = FVector2D::ZeroVector;
			int32 touchesLength = Corners[i]->Touches.Num();
			for (int j = 0; j < touchesLength; j++)
			{
				FMapCenter* center = Corners[i]->Touches[j];
				point += center->CenterData.Point;
			}
			point /= touchesLength;
			newCorners[i] = point;
		}
	}

	// Move the corners to the new locations.
	for (int i = 0; i < Corners.Num(); i++)
	{
		Corners[i]->CornerData.Point = newCorners[i];
	}

	// The edge midpoints were computed for the old corners and need
	// to be recomputed.
	for (int i = 0; i < Edges.Num(); i++)
	{
		if (Edges[i]->VoronoiEdge0 == NULL || Edges[i]->VoronoiEdge1 == NULL)
		{
			continue;
		}
		FVector2D edgeVertexOne = Edges[i]->VoronoiEdge0->CornerData.Point;
		FVector2D edgeVertexTwo = Edges[i]->VoronoiEdge1->CornerData.Point;
		Edges[i]->Midpoint = FVector2D(FMath::Lerp(edgeVertexOne.X, edgeVertexTwo.X, 0.5f), FMath::Lerp(edgeVertexOne.Y, edgeVertexTwo.Y, 0.5f));
	}
}

TArray<int32> UPolygonMap::FindLandCorners() const
{
	TArray<int32> landCorners;
	for (int i = 0; i < Corners.Num(); i++)
	{
		FMapCorner* corner = Corners[i];
		if (UMapDataHelper::IsOcean(corner->CornerData) || UMapDataHelper::IsCoast(corner->CornerData))
		{
			continue;
		}
		landCorners.Emplace(i);
	}
	return landCorners;
}

TArray<FMapData>& UPolygonMap::GetAllMapData()
{
	return CachedMapData;
}

void UPolygonMap::CompileMapData()
{
	CachedMapData.Empty();
	for (int i = 0; i < GetCenterNum(); i++)
	{
		CachedMapData.Add(GetCenter(i).CenterData);
	}
	for (int i = 0; i < GetCornerNum(); i++)
	{
		CachedMapData.Add(GetCorner(i).CornerData);
	}
}

FVector2D UPolygonMap::ConvertWorldPointToGraphSpace(const FVector WorldPoint, const FWorldSpaceMapData& WorldData, int32 MapSize)
{
	float xyScale = WorldData.XYScaleFactor / MapSize;

	FVector2D graphLocation = FVector2D::ZeroVector;
	graphLocation.X = WorldPoint.X / MapSize / xyScale;
	graphLocation.Y = WorldPoint.Y / MapSize / xyScale;

	return graphLocation;
}

FVector UPolygonMap::ConvertGraphPointToWorldSpace(const FMapData& MapData, const FWorldSpaceMapData& WorldData, int32 MapSize)
{
	float elevationOffset = WorldData.ElevationOffset;
	float xyScale = WorldData.XYScaleFactor / MapSize;
	float elevationScale = WorldData.ElevationScale;

	FVector worldLocation = FVector::ZeroVector;
	worldLocation.X = MapData.Point.X * MapSize * xyScale;
	worldLocation.Y = MapData.Point.Y * MapSize * xyScale;
	worldLocation.Z = (MapData.Elevation * elevationScale) + elevationOffset;

	return worldLocation;
}

FMapCenter& UPolygonMap::FindMapCenterForCoordinate(const FVector2D& Point) const
{
	if (Point.X > MaxPointLocation || Point.Y > MaxPointLocation || Point.X < MinPointLocation || Point.Y < MinPointLocation)
	{
		// Point out of bounds
		return *DefaultCenter;
	}
	for (int i = 0; i < Centers.Num(); i++)
	{
		FMapCenter* center = Centers[i];
		if (CenterContainsPoint(Point, *center))
		{
			return *center;
		}
	}
	return *DefaultCenter;
}

FMapCorner& UPolygonMap::FindMapCornerForCoordinate(const FVector2D& Point) const
{
	if (Point.X > MaxPointLocation || Point.Y > MaxPointLocation || Point.X < MinPointLocation || Point.Y < MinPointLocation)
	{
		// Point out of bounds
		return *DefaultCorner;
	}
	for (int i = 0; i < Corners.Num(); i++)
	{
		FMapCorner* corner = Corners[i];
		if (CornerContainsPoint(Point, *corner))
		{
			return *corner;
		}
	}
	return *DefaultCorner;
}

bool UPolygonMap::CenterContainsPoint(const FVector2D& Point, const FMapCenter& Center) const
{
	float minX = Center.CenterData.Point.X;
	float maxX = Center.CenterData.Point.X;
	float minY = Center.CenterData.Point.Y;
	float maxY = Center.CenterData.Point.Y;
	for (int i = 0; i < Center.Corners.Num(); i++)
	{
		// Iterate over the borders, creating a bounding box for the polygon
		FMapCorner corner = *Center.Corners[i];
		if(corner.Index < 0)
		{
			// Invalid corner
			continue;
		}
		maxX = FMath::Max(maxX, corner.CornerData.Point.X);
		minX = FMath::Min(minX, corner.CornerData.Point.X);
		maxY = FMath::Max(maxY, corner.CornerData.Point.Y);
		minY = FMath::Min(minY, corner.CornerData.Point.Y);
	}

	// Check to see if bounding box contains point
	if(Point.X < minX || Point.X > maxX || Point.Y < minY || Point.Y > maxY)
	{
		// Not even in the polygon's bounding box
		return false;
	}

	// Set the padding to be 1% of the polygon size
	float epsilon = ((maxX - minX) / 100.0f);
	FVector2D pointTarget = FVector2D(minX - epsilon, Point.Y);

	int32 intersections = 0;
	for (int i = 0; i < Center.Borders.Num(); i++)
	{
		FMapEdge* border = Center.Borders[i];
		if (border->Index < 0)
		{
			// Invalid index
			UE_LOG(LogWorldGen, Warning, TEXT("Invalid border! Center number %d, Edge index was %d."), Center.Index, Center.Borders[i]);
			continue;
		}
		if (SegementsIntersect(*border, Point, pointTarget))
		{
			intersections++;
			//UE_LOG(LogWorldGen, Log, TEXT("Checking if point (%f, %f)-(%f, %f) and point (%f, %f)-(%f, %f) intersect."), p1.X, p1.Y, q1.X, q1.Y, p2.X, p2.Y, q2.X, q2.Y);
		}
	}
	
	return (intersections & 1) == 1; // True if point is odd (inside of polygon)
}

bool UPolygonMap::CornerContainsPoint(const FVector2D& Point, const FMapCorner& Corner) const
{
	if (Corner.Touches.Num() != 3)
	{
		return false;
	}
	FVector2D p1 = Corner.Touches[0]->CenterData.Point;
	FVector2D p2 = Corner.Touches[1]->CenterData.Point;
	FVector2D p3 = Corner.Touches[2]->CenterData.Point;

	float y1 = p1.Y;
	float y2 = p2.Y;
	float y3 = p3.Y;

	float x1 = p1.X;
	float x2 = p2.X;
	float x3 = p3.X;

	// Calculate determinant
	float det = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3);
	if (det == 0.0f)
	{
		// Shouldn't happen, but okay
		return false;
	}

	float x = Point.X;
	float y = Point.Y;

	// https://stackoverflow.com/questions/36090269/finding-height-of-point-on-height-map-triangles
	float a = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / det;
	float b = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / det;
	float c = 1 - a - b;

	// p lies in T if and only if 0 <= a <= 1 and 0 <= b <= 1 and 0 <= c <= 1
	return 0 <= a && a <= 1 && 0 <= b && b <= 1 && 0 <= c && c <= 1;
}

// The main function that returns true if line segment 'p1q1'
// and 'p2q2' intersect.
bool UPolygonMap::SegementsIntersect(const FMapEdge& Edge, const FVector2D& StartPoint, const FVector2D& EndPoint) const
{
	FMapCorner* corner1 = Edge.VoronoiEdge0;
	FMapCorner* corner2 = Edge.VoronoiEdge1;
	if (corner1 == NULL || corner2 == NULL)
	{
		// Invalid corner
		UE_LOG(LogWorldGen, Warning, TEXT("Invalid corner! Edge number %d, Corner 1 index was %d, Corner 2 index was %d."), Edge.Index, Edge.VoronoiEdge0, Edge.VoronoiEdge1);
		return false;
	}
	const float v1x1 = corner1->CornerData.Point.X;
	const float v1y1 = corner1->CornerData.Point.Y;
	const float v1x2 = corner2->CornerData.Point.X;
	const float v1y2 = corner2->CornerData.Point.Y;
	const float v2x1 = StartPoint.X;
	const float v2y1 = StartPoint.Y;
	const float v2x2 = EndPoint.X;
	const float v2y2 = EndPoint.Y;

	float d1, d2;
	float a1, a2, b1, b2, c1, c2;

	// Convert the first vector to a line of infinite length
	// We want this line in a linear equation
	a1 = v1y2 - v1y1;
	b1 = v1x1 - v1x2;
	c1 = (v1x2 * v1y1) - (v1x1 * v1y2);

	// Every point (x, y) that solves the equation above is on the line.
	// Every point that does not solve it is not.
	// The equation will have a positive result if it is on one side of the line
	// and a negative result if it is on the other.
	// We insert our start and end points into the equation above.
	d1 = (a1 * v2x1) + (b1 * v2y1) + c1;
	d2 = (a1 * v2x2) + (b1 * v2y2) + c1;

	// If d1 and d2 have the same sign, they are both on the same side of our Edge.
	// In this case, no intersection is possible.
	// 0 is a special case, which is why < and > are used instead of <= and >=.
	if((d1 > 0 && d2 > 0) || (d1 < 0 && d2 < 0))
	{
		return false;
	}

	// If we have gotten to this point, the infinite lines intersect somewhere (i.e., they are not parallel).
	// However, we are looking for an intersection in a specific interval (the edge).
	// The edge is a subset of the infinite line, and it is possible that the intersection happened before or after that interval
	// In order to know for sure, we have to do the same test the other way around.
	// Starting with converting the second segment into an infinite line:
	a2 = v2y2 - v2y1;
	b2 = v2x1 - v2x2;
	c2 = (v2x2 * v2y1) - (v2x1 * v2y2);

	// Calculate d1 and d2 again, using the other segment
	d1 = (a2 * v1x1) + (b2 * v1y1) + c2;
	d2 = (a2 * v1x2) + (b2 * v1y2) + c2;

	// Again, if both have the same sign, no intersection is possible
	if ((d1 > 0 && d2 > 0) || (d1 < 0 && d2 < 0))
	{
		return false;
	}

	// Now we have 2 options.
	// Either they intersect at exactly one point, or they are colinear and intersect at any number of points.
	// We treat colinear lines as failing this test.
	if((a1 * b2) - (a2 * b1) < 0.001f)
	{
		// Points are colinear
		return false;
	}
	else
	{
		// Points intersect somewhere
		return true;
	}
}
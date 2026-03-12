#include "elasticPlate.h"

elasticPlate::elasticPlate(double m_YoungM_1, double m_density_1, double m_thickness_1, double m_Possion_1, 
		double m_YoungM_2, double m_density_2, double m_thickness_2, double m_Possion_2, 
		double m_dt)
{
	dt = m_dt;

	// material for first layer
	YoungM_1 = m_YoungM_1;
	density_1 = m_density_1;
	thickness_1 = m_thickness_1;
	Possion_1 = m_Possion_1;
	alpha_1 = YoungM_1 * Possion_1 / ( 1 - Possion_1 * Possion_1 );
	beta_1 = YoungM_1 / ( 2 * (1 + Possion_1) );

	// material for second layer
	YoungM_2 = m_YoungM_2;
	density_2 = m_density_2;
	thickness_2 = m_thickness_2;
	Possion_2 = m_Possion_2;
	alpha_2 = YoungM_2 * Possion_2 / ( 1 - Possion_2 * Possion_2 );
	beta_2 = YoungM_2 / ( 2 * (1 + Possion_2) );

	setupGeometry();

	ndof = 3 * nv;
	x = VectorXd::Zero(ndof);
	x0 = VectorXd::Zero(ndof);
	u = VectorXd::Zero(ndof);

	for (int i = 0; i < nv; i++)
	{
		x(3 * i + 0) = v_nodes[i](0);
		x(3 * i + 1) = v_nodes[i](1);
		x(3 * i + 2) = v_nodes[i](2);
	}
	x0 = x;

	initializeTriangular();

	setupMass();

	//set up constraint map
	isConstrained = new int[ndof];
    for (int i=0; i < ndof; i++)
    {
		isConstrained[i] = 0;
    }
}



elasticPlate::~elasticPlate()
{
	delete isConstrained;
	delete unconstrainedMap;
	delete fullToUnconsMap;
}

void elasticPlate::setup()
{
	ncons = 0;
    for (int i=0; i < ndof; i++)
    {
		if (isConstrained[i] > 0)
		{
			ncons++;
		}
	}
	uncons = ndof - ncons;

	unconstrainedMap = new int[uncons]; // maps xUncons to x
	fullToUnconsMap = new int[ndof];
	setupMap();
}

void elasticPlate::setupMap()
{
	int c = 0;
	for (int i=0; i < ndof; i++)
	{
		if (isConstrained[i] == 0)
		{
			unconstrainedMap[c] = i;
			fullToUnconsMap[i] = c;
			c++;
		}
	}
}

void elasticPlate::setupMass()
{
	massArray = VectorXd::Zero(ndof);

	for (int i = 0; i < triangularNum; i++)
	{
		int index1 = v_triangularElement[i].nv_1;
		int index2 = v_triangularElement[i].nv_2;
		int index3 = v_triangularElement[i].nv_3;

		double deltaMass = v_triangularElement[i].area * ( thickness_1 * density_1 + thickness_2 * density_2 ) / 3;

		massArray(3 * index1 + 0) = massArray(3 * index1 + 0) + deltaMass;
		massArray(3 * index1 + 1) = massArray(3 * index1 + 1) + deltaMass;
		massArray(3 * index1 + 2) = massArray(3 * index1 + 2) + deltaMass;
		
		massArray(3 * index2 + 0) = massArray(3 * index2 + 0) + deltaMass;
		massArray(3 * index2 + 1) = massArray(3 * index2 + 1) + deltaMass;
		massArray(3 * index2 + 2) = massArray(3 * index2 + 2) + deltaMass;

		massArray(3 * index3 + 0) = massArray(3 * index3 + 0) + deltaMass;
		massArray(3 * index3 + 1) = massArray(3 * index3 + 1) + deltaMass;
		massArray(3 * index3 + 2) = massArray(3 * index3 + 2) + deltaMass;
	}
}

int elasticPlate::getIfConstrained(int k)
{
	return isConstrained[k];
}

void elasticPlate::setVertexBoundaryCondition(Vector3d position, int k)
{
	isConstrained[3 * k + 0] = 1;
	isConstrained[3 * k + 1] = 1;
	isConstrained[3 * k + 2] = 1;

	// Store in the constrained dof vector
	x(3 * k + 0) = position(0);
	x(3 * k + 1) = position(1);
	x(3 * k + 2) = position(2);
}

void elasticPlate::setConstraint(double position, int k)
{
	isConstrained[k] = 1;

	// Store in the constrained dof vector
	x(k) = position;
}

void elasticPlate::setOneVertexBoundaryCondition(double position, int i, int k)
{
	isConstrained[3 * i + k] = 1;

	// Store in the constrained dof vector
	x(3 * i + k) = position;
}

Vector3d elasticPlate::getVertex(int i)
{
	Vector3d xCurrent;

	xCurrent(0) = x(3 * i + 0);
	xCurrent(1) = x(3 * i + 1);
	xCurrent(2) = x(3 * i + 2);

	return xCurrent;
}

Vector3d elasticPlate::getVertexOld(int i)
{
	Vector3d xCurrent;

	xCurrent(0) = x0(3 * i + 0);
	xCurrent(1) = x0(3 * i + 1);
	xCurrent(2) = x0(3 * i + 2);

	return xCurrent;
}

Vector3d elasticPlate::getVelocity(int i)
{
	Vector3d uCurrent;

	uCurrent(0) = ( x(3 * i + 0) - x0(3 * i + 0) ) / dt;
	uCurrent(1) = ( x(3 * i + 1) - x0(3 * i + 1) ) / dt;
	uCurrent(2) = ( x(3 * i + 2) - x0(3 * i + 2) ) / dt;

	uCurrent(0) = u(3 * i + 0);
	uCurrent(1) = u(3 * i + 1);
	uCurrent(2) = u(3 * i + 2);

	return uCurrent;
}

void elasticPlate::initializeTriangular()
{
	triangularNum = 0;
	v_triangularElement.clear();

	for (int i = 0; i < triangular.size(); i++)
	{
		Vector3i triangularCurrent = triangular[i];

		triangularElement m_triangularElement;

		m_triangularElement.nv_1 = triangularCurrent(0);
		m_triangularElement.nv_2 = triangularCurrent(1);
		m_triangularElement.nv_3 = triangularCurrent(2);

		m_triangularElement.x_1 = getVertex(m_triangularElement.nv_1);
		m_triangularElement.x_2 = getVertex(m_triangularElement.nv_2);
		m_triangularElement.x_3 = getVertex(m_triangularElement.nv_3);

		// m_triangularElement.alpha = alpha;
		// m_triangularElement.beta  = beta;
		// m_triangularElement.thickness = thickness;

		Vector3d x_1 = getVertex(m_triangularElement.nv_1);
		Vector3d x_2 = getVertex(m_triangularElement.nv_2);
		Vector3d x_3 = getVertex(m_triangularElement.nv_3);

		Vector3d e_1 = x_2 - x_1;
		Vector3d e_2 = x_3 - x_1;
		m_triangularElement.area = 0.5 * ( e_1.cross(e_2) ).norm();

		m_triangularElement.abar_1 = getFFF(i, nullptr, nullptr);
		
		m_triangularElement.abarinv_1 = m_triangularElement.abar_1.inverse();

		m_triangularElement.abar_2 = m_triangularElement.abar_1;
		m_triangularElement.abarinv_2 = m_triangularElement.abarinv_1;

		// double detA = m_triangularElement.abar_1.determinant();
		

		// b bar

		
		Matrix2d B = getSFF(i, nullptr, nullptr);
		
		m_triangularElement.bbar_1 = B;

        m_triangularElement.bbar_2 = m_triangularElement.bbar_1;

		m_triangularElement.arrayIndex = VectorXi::Zero(9);

		m_triangularElement.arrayIndex(0) = 3 * m_triangularElement.nv_1 + 0;
		m_triangularElement.arrayIndex(1) = 3 * m_triangularElement.nv_1 + 1;
		m_triangularElement.arrayIndex(2) = 3 * m_triangularElement.nv_1 + 2;

		m_triangularElement.arrayIndex(3) = 3 * m_triangularElement.nv_2 + 0;
		m_triangularElement.arrayIndex(4) = 3 * m_triangularElement.nv_2 + 1;
		m_triangularElement.arrayIndex(5) = 3 * m_triangularElement.nv_2 + 2;

		m_triangularElement.arrayIndex(6) = 3 * m_triangularElement.nv_3 + 0;
		m_triangularElement.arrayIndex(7) = 3 * m_triangularElement.nv_3 + 1;
		m_triangularElement.arrayIndex(8) = 3 * m_triangularElement.nv_3 + 2;

		v_triangularElement.push_back(m_triangularElement);

		triangularNum = triangularNum + 1;		
	}
}

void elasticPlate::prepareForIteration()
{
	;
}

void elasticPlate::updateTimeStep()
{
	prepareForIteration();

	// compute velocity
	u = (x - x0) / dt;

	// update x
	x0 = x;
}

void elasticPlate::updateGuess()
{
	for (int c=0; c < uncons; c++)
	{
		x[unconstrainedMap[c]] = x[unconstrainedMap[c]] + u[unconstrainedMap[c]] * dt;
	}
}

void elasticPlate::updateNewtonMethod(VectorXd m_motion)
{
	for (int c=0; c < uncons; c++)
	{
		x[unconstrainedMap[c]] -= m_motion[c];
	}
}

void elasticPlate::setupGeometry()
{
	nv = 0;
	triangularNum = 0;

	ifstream inFile1;
	inFile1.open("inputdata/nodesInput.txt");
	v_nodes.clear();
	double a, b, c;
	while(inFile1 >> a >> b >> c)
	{
		Vector3d xCurrent;

		xCurrent(0) = a;
		xCurrent(1) = b;
		xCurrent(2) = c;

		v_nodes.push_back(xCurrent);
	}
	nv = v_nodes.size();
	inFile1.close();

	ifstream inFile2;
	inFile2.open("inputdata/triangleInput.txt");
	v_triangularElement.clear();
	int aa, bb, cc;
	while(inFile2 >> aa >> bb >> cc)
	{
		aa = aa;
		bb = bb;
		cc = cc;

		Vector3i triangularCurrent;
		triangularCurrent(0) = aa;
		triangularCurrent(1) = bb;
		triangularCurrent(2) = cc;
		triangular.push_back(triangularCurrent);
	}
	triangularNum = triangular.size();
	inFile2.close();

    F = MatrixXi::Zero(triangularNum, 3);

    for(int i = 0; i < triangular.size(); i++)
    {
        Vector3i triangularCurrent = triangular[i];

        F.row(i) = triangularCurrent;
    }

    buildMeshConnectivity(F);

    // std::cout << "nfaces = " << mesh.nfaces << "\n";
    // std::cout << "nedges = " << mesh.nedges << "\n";
    // std::cout << "EV:\n" << mesh.EV << "\n";
    // std::cout << "EF:\n" << mesh.EF << "\n";
    // std::cout << "EOpp:\n" << mesh.EOpp << "\n";
    // std::cout << "FE:\n" << mesh.FE << "\n";
    // std::cout << "FEorient:\n" << mesh.FEorient << "\n";
    // std::cout << "Fbendvert:\n" << mesh.Fbendvert << "\n";
    // std::cout << "Fbenddof:\n" << mesh.Fbenddof << "\n";
    // std::cout << "isVirtualNormal:\n" << mesh.isVirtualNormal << "\n";

    totalEdge = mesh.nedges;
}

void elasticPlate::buildMeshConnectivity(MatrixXi F) 
{
		mesh.F = F;

    	std::map<std::pair<int, int>, Eigen::Vector2i > edgeFaces;

        int nfaces = (int)F.rows();
        mesh.nfaces = nfaces;

        
        for (int i = 0; i < nfaces; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                int v0 = F(i, (j + 1) % 3);
                int v1 = F(i, (j + 2) % 3);
                int idx = 0;
                if (v0 > v1)
                {
                    std::swap(v0, v1);
                    idx = 1;
                }

                std::pair<int, int> p(v0, v1);
                auto it = edgeFaces.find(p);
                if (it == edgeFaces.end())
                {
                    edgeFaces[p][idx] = i;
                    edgeFaces[p][1 - idx] = -1;
                }
                else
                {
                    edgeFaces[p][idx] = i;
                }
            }
        }

        

        int nedges = (int)edgeFaces.size();
        mesh.nedges = nedges;

        mesh.FE.resize(nfaces, 3);
        mesh.FEorient.resize(nfaces, 3);
        mesh.EV.resize(nedges, 2);
        mesh.EF.resize(nedges, 2);
        mesh.EOpp.resize(nedges, 2);
        std::map<std::pair<int, int>, int> edgeIndices;

        int idx = 0;
        for (auto it : edgeFaces)
        {
            edgeIndices[it.first] = idx;
            mesh.EV(idx, 0) = it.first.first;
            mesh.EV(idx, 1) = it.first.second;

            if (it.second[0] < it.second[1] )
            {
            	mesh.EF(idx, 0) = it.second[0];
            	mesh.EF(idx, 1) = it.second[1];
            }
            else
            {
            	mesh.EF(idx, 0) = it.second[1];
            	mesh.EF(idx, 1) = it.second[0];
            }

            if (mesh.EF(idx, 0) < 0)
            {
            	int temp = mesh.EF(idx, 0);
            	mesh.EF(idx, 0) = mesh.EF(idx, 1);
            	mesh.EF(idx, 1) = temp;
            }
       
            idx++;
        }

        for (int i = 0; i < nfaces; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                int v0 = F(i, (j + 1) % 3);
                int v1 = F(i, (j + 2) % 3);
                if (v0 > v1) std::swap(v0, v1);
                mesh.FE(i, j) = edgeIndices[std::pair<int, int>(v0, v1)];
            }
        }

        for (int i = 0; i < nedges; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                mesh.EOpp(i, j) = oppositeVertex(i, j);
            }
        }

        for (int i = 0; i < nfaces; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                int edgeID = mesh.FE(i, j);
        		Vector2i edgeVerts = mesh.EV.row(edgeID);
        
        		int v1 = mesh.F(i, (j+1) % 3);
        		int v2 = mesh.F(i, (j+2) % 3);

        		if ( v1 == edgeVerts(0) && v2 == edgeVerts(1) )
        		{
        			mesh.FEorient(i, j) = 0;
        		}
        		else
        		{
        			mesh.FEorient(i, j) = 1;
        		}
            }
        }

        mesh.Fbendvert = MatrixXi::Constant(nfaces, 6, -1);
    	mesh.Fbenddof  = MatrixXi::Constant(nfaces, 18, -1);

    	for (int i = 0; i < nfaces; ++i) 
    	{
        	std::vector<int> verts = {F(i, 0), F(i, 1), F(i, 2)};

        	for (int j = 0; j < 3; ++j) 
        	{
            	mesh.Fbendvert(i, j) = verts[j];
            	mesh.Fbenddof(i, j * 3 + 0) = 3 * verts[j] + 0;
            	mesh.Fbenddof(i, j * 3 + 1) = 3 * verts[j] + 1;
            	mesh.Fbenddof(i, j * 3 + 2) = 3 * verts[j] + 2;
        	}

        	for (int j = 0; j < 3; ++j) 
        	{
            	int edgeID = mesh.FE(i, j);
            	std::vector<int> eOppVerts = {mesh.EOpp(edgeID, 0), mesh.EOpp(edgeID, 1)};
            	// Remove any that are in face verts
            	for (int v : verts) 
            	{
                	eOppVerts.erase(std::remove(eOppVerts.begin(), eOppVerts.end(), v), eOppVerts.end());
           	 	}
            	if (!eOppVerts.empty()) 
            	{
                	int opp = eOppVerts[0];
                	mesh.Fbendvert(i, 3 + j) = opp;

                	mesh.Fbenddof(i, 9 + j * 3 + 0) = 3 * opp + 0;
               	 	mesh.Fbenddof(i, 9 + j * 3 + 1) = 3 * opp + 1;
                	mesh.Fbenddof(i, 9 + j * 3 + 2) = 3 * opp + 2;
            }
        }
    }

    mesh.isVirtualNormal = VectorXi::Zero(mesh.nedges);  
}

Vector3d elasticPlate::faceNormal(const Vector3d& q0, const Vector3d& q1, const Vector3d& q2)
{
	Vector3d n = (q1 - q0).cross(q2 - q0);

	return n;
}

int elasticPlate::edgeVertex(int edge, int vertidx) 
{ 
	return mesh.EV(edge, vertidx); 
}

int elasticPlate::edgeFace(int edge, int faceidx)
{ 
	return mesh.EF(edge, faceidx); 
}

int elasticPlate::edgeOppositeVertex(int edge, int faceidx) 
{ 
	return mesh.EOpp(edge, faceidx); 
}

int elasticPlate::faceVertex(int face, int vertidx)
{ 
	return mesh.F(face, vertidx); 
}

int elasticPlate::faceEdge(int face, int vertidx)
{ 
	return mesh.FE(face, vertidx);
} 

int elasticPlate::faceEdgeOrientation(int face, int vertidx)
{
	return mesh.FEorient(face, vertidx);
} 

int elasticPlate::oppositeVertexIndex(int edge, int faceidx)
{
    int face = edgeFace(edge, faceidx);

    if (face == -1)
    {
        return -1;
    }

    for (int j = 0; j < 3; j++)
    {
        if (F(face, j) != edgeVertex(edge, 0) && F(face, j) != edgeVertex(edge, 1))
        return j;
    }

    return -1;
}

int elasticPlate::oppositeVertex(int edge, int faceidx)
{
    int face = edgeFace(edge, faceidx);
    int idx = oppositeVertexIndex(edge, faceidx);

    if (idx == -1)
    {
        return -1;
    }

    return F(face, idx);
}

int elasticPlate::vertexOppositeFaceEdge(int face, int vertidx)
{
    int edge = faceEdge(face, vertidx);
    int edgeorient = faceEdgeOrientation(face, vertidx);

    return edgeOppositeVertex(edge, 1 - edgeorient);
}

Matrix2d elasticPlate::getFFF(int idx , 
	Matrix<double, 4, 9>* derivative,
    std::vector<Matrix<double, 9, 9> >* hessian)
{
	
	Vector3i cfaceidx = mesh.F.row(idx);
	Vector3d q0 = getVertex(cfaceidx(0));
	Vector3d q1 = getVertex(cfaceidx(1));
	Vector3d q2 = getVertex(cfaceidx(2));
	Matrix2d A = Geometry::firstFundamentalForm( q0 , q1 , q2 , derivative,  hessian );

	return A;
}


Matrix2d elasticPlate::getSFF(int idx , 
	Matrix<double, 4, 18>* derivative,
    std::vector<Matrix<double, 18, 18> >* hessian)
{
	// build sFF information for the face
	Vector3d x_1 = getVertex(mesh.F(idx,0));
	Vector3d x_2 = getVertex(mesh.F(idx,1));
	Vector3d x_3 = getVertex(mesh.F(idx,2));

	Geometry::sFFinformation sFFinf(x_1, x_2, x_3);
	for (int v = 0; v < 3; v++)
	{
		
		int adjnode = mesh.Fbendvert(idx, v + 3);
		int edge = mesh.FE(idx, v);		// The adjnode and edge are corresponding, adjnode is the opposite vertex of the edge.
		if ( adjnode != -1)
		{ 
			Vector3d position =  getVertex(adjnode);
			sFFinf.setupAdjVertex(v, position);
		}
		else if (mesh.isVirtualNormal(edge) != 0) 	//setup virtual normal if need (for future)
		{
			Vector3d normal = {0, 0, 0};
			sFFinf.setupVirtualNormal(v, normal);
		}
		
	}
	Matrix2d B = Geometry::secondFundamentalForm(sFFinf, derivative, hessian);

	return B;
}


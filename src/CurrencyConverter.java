import java.util.Map;
import java.util.List;
import java.util.Set;
import java.util.HashSet;
import java.util.HashMap;
import java.util.Arrays;
import java.util.ArrayList;

//over all TC: O((V+E)*log(V)) bc of dijkstra
//overall SC: O(V+E) bc of dijkstra
public class CurrencyConverter {


    /**
     * Finds the path that leads to the max amount of money.
     * TC: O((V+E) * log(V)) where V is vertices and E is edges bc you need to modifiy edges in edgeList
     * time complexity of Dij algo is O((V+E)*log(V)) bc you update a distance each time, at TC of findingConversionPath
     * SC: O(V+E) where V is number of Vertices and E is number of edges bc you need to store number of modified edges
     * and dij algo uses O(V) space for converted distances map, and SC of FindingConversionPath
     *
     * @param edgeList
     * @param fromCurrency
     * @param toCurrency
     */
    public static void dijkstra(List<List<Object>> edgeList, Object fromCurrency, Object toCurrency){

        //modifiy edgeList to use neg logs of exchange rates as weights
        List<List<Object>> modifiedEdgeList = new ArrayList<>();
        for(List<Object> edge : edgeList){
            double weight = -Math.log((double) edge.get(2)); //Negative logs of exchange rates
            modifiedEdgeList.add(Arrays.asList(edge.get(0), edge.get(1), weight));
        }

        //using dijkstra's algo with modified edge list
        Map<Object, Double> convertedDistances = runDijkstra(modifiedEdgeList, fromCurrency);
        //finding the converstion path the leads to the most money using dfs
        List<Object> path = findConversionPath(modifiedEdgeList, convertedDistances, fromCurrency, toCurrency);

        //print out result
        System.out.println("Path to most money: " + path);
    }

    /**
     * Uses a modification of Dijkstra's Algorithm to
     * convert edges into negative logs.
     * TC: O((V-1)*E) where V is number of vertices and E is the number
     * of edges because edges are iterated over once O(E), and vertices are iterated over
     * once to create the adj matrix m O(V), and initalizing v's as inf takes O(V).
     * Running through Dij algo: outer loop O(V-1), inner its. over edgList O(E).
     * Converting distances take O(V)
     *
     * SC: O(V) where V is number of vertices because vertices set stores vertices from edgeList
     * VertexIndexMap stores vertices and indices & distance stores distance from source to vertices
     * ConvertedIDstances is also O(V) bc it converts di for eac v
     * @param edgeList
     * @param source
     * @return
     */
    //Object is a generic type that can rep any java object
    public static Map<Object, Double> runDijkstra(List<List<Object>> edgeList, Object source){

        //getting vertices from edgeList, iterates ONCE (set contains no dups)
        Set<Object> vertices = new HashSet<>();
        for(List<Object> edge : edgeList){
            vertices.add(edge.get(0)); //index of first ele in list (starting vertex)
            vertices.add(edge.get(1)); //index of second ele in list (ending vertex)
        }

        //creating an adj matrix (each vertex is mapped to an index)
        Map<Object, Integer> vertexIndexMap = new HashMap<>(); //stores mapping between vertices & corresponding indices
        int index = 0; //starting at 0 to assign indices to each vertex
        for(Object vertex : vertices){
            vertexIndexMap.put(vertex, index++);
        }

        //initializing array to store distances from source to all other vertices
        int v = vertices.size();
        double[] distance = new double[v]; //stores the distances
        Arrays.fill(distance, Double.POSITIVE_INFINITY); //setting them to infinity
        distance[vertexIndexMap.get(source)] = 0; //sets start vertex at 0

        //main part of dijkstra's algorithm
        for(int i = 0; i < v -1; i++){ //iterate over vertices (worst case: visiting all vertices)

            for(List<Object> edge : edgeList){ //iterates over the edges

                int fromIndex = vertexIndexMap.get(edge.get(0)); //gets start vertex of current edge
                int toIndex = vertexIndexMap.get(edge.get(1)); //ending vertex of current edge
                double weight = (double) edge.get(2); //gets weight of current edge

                //checks if the distance from source vertex to ending vertex
                //is less than the current shortest distance stored in distance[toIndex]
                if(distance[fromIndex] + weight < distance[toIndex]){

                    /*
                    if true, then update the shortest diatnce from the source vertex
                    to the ending vertex of the current edge to the sum of the distances
                    from the source vertex to the strting vertex of the current edge
                    and the weight of the current edge and update the shortest path if a
                    shorter path is found
                     */
                    distance[toIndex] = distance[fromIndex] + weight;
                }
            }
        }

        Map<Object, Double> convertedDistances = new HashMap<>(); //stores convertedDistances for each vertex
        //loops over vertexIndexMap and which has their corresponding vertices and indeices
        for(Map.Entry<Object, Integer> entry : vertexIndexMap.entrySet()){

            /*
            calculates the converted distances for the current vertex
            gets the current vertex's distance and converts it into a neg log
             */
            convertedDistances.put(entry.getKey(), Math.exp(-distance[entry.getValue()]));
        }

        return convertedDistances;
    }

    /**
     * TC: O(V+E) where V is the number of vertices and E is the number of edges
     *  in the WC be dfs msut visited each vertex and edge once.
     *  SC: O(V+E) where V is number of vertices and E is number of edges bc space
     *  allocated for ArrayList, HashSet, and recursive stack depends on V and E.
     * @param edgeList
     * @param convertedDistances
     * @param source
     * @param destination
     * @return
     */
    private static List<Object> findConversionPath(List<List<Object>> edgeList, Map<Object, Double> convertedDistances,
                                                   Object source, Object destination){
        List<Object> path = new ArrayList<>(); //stores the conversion path
        Set<Object> visited = new HashSet<>(); //keep track of visited vertices
        dfs(edgeList, convertedDistances, source, destination, visited, path);
        return path; //returns the path list which contains conversion path from source currency to destination currency
    }

    /**
     * Uses dfs recursively to find a conversion path between
     * a given source and destination v's in a graph.
     * TC: O(E) (in WC) where E is number of edges because each edge is travered once using dfs.
     * In WC if graph is represented as adj matrix and all vertices are connected to each other
     * than TC can be O(V+E)
     * SC: O(V) where V is number of vertices becasue of the recursive call stack and data structures.
     * - recursive call stack depends on max depth of DFS
     * - auxillary data structures visited set & path list are equal to number of vertices in graph
     * @param edgeList
     * @param convertedDistances
     * @param current
     * @param destination
     * @param visited keep track of visited vertices
     * @param path stores the conversion path
     * @return the resulting conversion path
     */
    private static boolean dfs(List<List<Object>> edgeList, Map<Object, Double> convertedDistances,
                               Object current, Object destination, Set<Object> visited, List<Object> path){

        //checks if current vertex is desination vertex
        if(current.equals(destination)){
            //if true, adds current to vertex path list and returns true indicating destination is reached
            path.add(current);
            return true;
        }

        //otherwise mark current vertex as visited
        visited.add(current);
        //iterate through edgeList to find edges starting from current vertex
        for(List<Object> edge : edgeList){
            //for each edge, if starting vertex is the current & ending vertex isn't visited
            //recursviely invokes dfs with ending vertex
            if(edge.get(0).equals(current) && !visited.contains(edge.get(1))){
                //if recursive call returns true, (meaning destination is reached)
                //it adds the current vertex to begining of path and returns true
                if(dfs(edgeList, convertedDistances, edge.get(1), destination, visited, path)){
                    path.add(0, current);
                    return true;
                }
            }
        }
        return false; //if no path is found from current to destination using dfs it returns false
    }

    public static void main(String[] args){

        //test 1
        List<List<Object>> edgeList1 = new ArrayList<>();
        edgeList1.add(Arrays.asList("USD", "JPY", 110.0));
        edgeList1.add(Arrays.asList("USD", "AUD", 1.45));
        edgeList1.add(Arrays.asList("JPY", "COP", 1.3));
        System.out.println("Test Case 1:");
        dijkstra(edgeList1, "USD", "COP");

        System.out.println();
        List<List<Object>> edgeList2 = new ArrayList<>();
        edgeList2.add(Arrays.asList("EUR", "USD", 1.2));
        edgeList2.add(Arrays.asList("USD", "JPY", 110.0));
        edgeList2.add(Arrays.asList("JPY", "GBP", 0.007));
        edgeList2.add(Arrays.asList("EUR", "GBP", 0.9));
        System.out.println("Test Case 2:");
        dijkstra(edgeList2, "EUR", "JPY");
    }
}

package ComponentModels;

import java.awt.Color;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import Distributions.NormalDistribution;
import Graphics.GenGraphics;
import Mathematics.MatrixOperations;

public class TestICA {

	
	public static void test1() {
		
		//Plot: obs, ica, pca		
		String plotData = "pca";
		
		ArrayList<double [][]> org_signals = get_orginial_signals(false);
		HashMap<String, ArrayList<double [][]>> obs_signalList = get_observed_signals(org_signals);
		ArrayList<double [][]> obs_signals = obs_signalList.get("list");
		double [][] obs_signal_X = obs_signalList.get("array").get(0);
		
		ICA ica = new ICA(obs_signal_X, 4, "cosh", true, true);
		ica.do_ICA();
		double [][] icaReconstruction = ica.calc_uncentered_and_unwhitened_rotated_X();
		double [][] icaFactors = ica.get_factors();
		
		int n_variables = icaReconstruction[0].length;
		ArrayList<double [][]> icaReconsList = new ArrayList<double [][]>();
		ArrayList<double [][]> icaFactorList = new ArrayList<double [][]>();
		for(int i=0; i<n_variables; i++) {
			icaReconsList.add(MatrixOperations.get_column_vec_from_matrix(icaReconstruction,i));
			icaFactorList.add(MatrixOperations.get_column_vec_from_matrix(icaFactors,i));
		}
	
		PCA pca = new PCA(obs_signal_X, 4, false, false);
		pca.do_PCA();
		double [][] pcaReconstruction = pca.get_rotated_input();
		double [][] pcaFactors = pca.get_factors();

		n_variables = pcaReconstruction[0].length;
		ArrayList<double [][]> pcaReconsList = new ArrayList<double [][]>();
		ArrayList<double [][]> pcaFactorList = new ArrayList<double [][]>();
		for(int i=0; i<n_variables; i++) {
			pcaReconsList.add(MatrixOperations.get_column_vec_from_matrix(pcaReconstruction,i));
			pcaFactorList.add(MatrixOperations.get_column_vec_from_matrix(pcaFactors,i));
		}
	
		if(plotData.contentEquals("ica")) {
			plot_ica(icaFactorList,icaReconsList);
		}
		if(plotData.contentEquals("pca")) {
			plot_pca(pcaFactorList,pcaReconsList);
		}
		if(plotData.contentEquals("obs")) {
			plotSignals(org_signals, obs_signals, null);
		}
	}
	
	
	@SuppressWarnings("static-access")
	public static void test2() {
		
		int n=100;
		int max = 10;
		int min = -10;
		double diff = (max-min)/((double)(n));
		
		ArrayList<double [][]> org_signals = get_orginial_signals(false);
		HashMap<String, ArrayList<double [][]>> obs_signalList = get_observed_signals(org_signals);
		double [][] obs_signal_X = obs_signalList.get("array").get(0);
		
		double [][] x = new double [n][1];
		double [][] dist = new double [n][1];
		
		ICA ica = new ICA(obs_signal_X, 4, true, true);
		for(int i=0; i<n; i++) {
			x[i][0] = min+diff*i;
			dist[i][0] = ica.g_superGaussian(x[i][0], 0);
		}
		
		plotDistribution4ICA(x, dist, "Sub Gaussian");
		
	}
	

	public static HashMap<String, ArrayList<double [][]>> get_observed_signals(ArrayList<double [][]> org_signals) {
		
		HashMap<String, ArrayList<double [][]>> obs_signalList = new HashMap<String, ArrayList<double [][]>>();
		
		int n_vars = org_signals.size();
		int sampleLength = org_signals.get(0).length;
		
		double [][] mixing_matrix = new double [4][4];
		for(int i=0; i<n_vars; i++) {
			for(int j=0; j<n_vars; j++) {
				if(i==j) {
					mixing_matrix[i][j] = 1.5;
				}else {			
					if(i==2) {
						//Weighting of the white noise term
						mixing_matrix[i][j] = 1.2;
					}else {
						mixing_matrix[i][j] = 0.8;
					}														
				}				
			}
		}
		
		double [][] obs_X = new double [sampleLength][n_vars];
		for(int i=0; i<n_vars; i++) {
			for(int j=0; j<sampleLength; j++) {
				obs_X[j][i] = org_signals.get(i)[j][0];
			}
		}
		
		obs_X = MatrixOperations.multiplication(obs_X, mixing_matrix);
		
		ArrayList<double [][]> obs_signals = new ArrayList<double [][]>();
		for(int i=0; i<n_vars; i++) {
			obs_signals.add(MatrixOperations.get_column_vec_from_matrix(obs_X, i));
		}
		
		ArrayList<double [][]> arrayList = new ArrayList<double [][]>();
		arrayList.add(obs_X);
		
		obs_signalList.put("list", obs_signals);
		obs_signalList.put("array", arrayList);		
		
		return obs_signalList;
	}
	
	
	public static ArrayList<double [][]> get_orginial_signals(boolean createNewGaussianProcess) {
		
		int signalLength = 501;
		
		ArrayList<double [][]> org_signals = new ArrayList<double [][]>();
		
		double [][] sawtoothWave = sawtoothWave(signalLength, 25, 4.0, -2.0);		
		org_signals.add(sawtoothWave);
		
		double [][] tangWave = tangensWave(signalLength, 20, 4.0, 0.0);
		org_signals.add(tangWave);
		
		double [][] whiteNoise = new double [signalLength][1];
		if(createNewGaussianProcess == true) {
			whiteNoise = whiteNoiseProcess(signalLength);
		}else {
			whiteNoise = get_frozenWhiteNoiseProcess();
		}
		org_signals.add(whiteNoise);
		
		double [][] sinusWave = sinusWave(signalLength, 40, 1.0, 0.0);
		org_signals.add(sinusWave);
		
		return org_signals;
	}
	
	
	public static double [][] whiteNoiseProcess(int signalLength) {
		
		NormalDistribution nDist = new NormalDistribution(0.0, 2.0);
		double [][] whiteNoise = nDist.sample(signalLength);
		getWhiteNoiseAsString(whiteNoise);
		return whiteNoise;
	}
	
	
	public static void getWhiteNoiseAsString(double [][] whiteNoise){
		
		int n =  whiteNoise.length;
		String noise = ""; 
		for(int i=0; i<n; i++) {
			noise += whiteNoise[i][0];
			if(i<(n-1)) {
				noise += ", ";
			}
		}
		System.out.println(noise);
	}
	
	
	public static double [][] get_frozenWhiteNoiseProcess() {
		
		double [] whiteNoise = {1.5848219477918755, 3.1308638403164726, 0.397038237955182, -1.0409087108404884, 
								0.7119408863553193, 0.06969058470200455, 2.8496213193034086, 1.518680482315091, 
								0.5745327740575117, -0.23363998930992655, 0.02997072514917987, -0.8194952732540463, 
								1.4459028318913576, -0.3546620395099807, -2.629951076242764, 0.6029182779365321, 
								-1.5753128629395579, -1.4161841358147247, -1.0302673076986752, 1.4890037242359897, 
								1.7664385352576464, 1.3686539877904356, -1.6285625634707266, -1.8554844772201173, 
								0.6532016772606734, -1.3162581655148602, 0.7129559612243346, -0.9373398125439591, 
								0.5746990002587122, 1.7758930423851493, 0.5392022469879423, -0.6242924887532644, 
								-0.6410869649383971, 0.6268244964373935, -1.9113660620801738, -0.25633922902888845, 
								1.9861219750539822, 0.39004282062306683, -1.6883541091923648, 0.17592536575870943, 
								-1.0115190758179027, -0.23988455956368526, -0.5894643916538324, -2.9809970546508566, 
								0.23009908557822825, -0.6772150466235198, 1.5787676797544672, -0.413543145281927, 
								-1.0430354596880784, 1.7983692658833812, -1.0084610838973287, -1.6907985353739639, 
								-1.626456860167883, 0.248523115343813, -1.1791974539535148, -1.429021195463105, 
								1.7308187584670522, 1.01339924765421, -1.9379489497265863, -0.3945719217931045, 
								0.19117197955237583, 1.735304583140889, -3.5420752196630487, -1.205583620174505, 
								0.44549837481376425, -0.3199787213678959, 0.19053138884452198, -0.4091797418905186, 
								1.2808867131184924, 0.009858330148418611, -1.2352176638310797, 0.4008744146113583, 
								1.0549509099149674, -1.5734711315456027, -2.094200085030608, -0.1119696171669894, 
								0.4388648750137345, -0.3509968277818327, 0.7550363148107998, -1.6582802252655415, 
								0.7977500644283698, 0.7401335943278629, 0.12055891298026003, -0.0737553688935073, 
								0.5289569759647246, -0.03918815192996235, 3.0458610517092737, 0.9267213076576549, 
								1.4276069181685564, 2.1738373973591973, 0.4307346277057418, 2.9331013920224684, 
								-1.8315694074372, 0.1457477585615914, -0.12904271080089538, 0.9951473296832313, 
								-2.7948119455225107, -0.5843069375883455, -1.5630658765023666, 1.6800035143528205, 
								1.3018428528496082, 0.9225224567600692, -0.07390457398346902, -0.37219054801959683, 
								0.7669974761015065, 0.6686123916154978, 0.2635271961112167, -1.6618790973122584, 
								-1.7441070174190072, 0.2794818379572095, 1.4964925107245042, -0.8539209904051402, 
								-0.8495220604667941, 0.7810301574696779, -1.5218684182708015, -1.4512747085986046, 
								0.7542224367535468, -0.3041300882959847, -1.377636060803164, -1.7779061632082132, 
								0.4193105631828722, -0.5180724644172008, 1.3021913794120676, -0.7679178663540153, 
								-2.2928609700390616, -1.7355181426476354, -1.3929701884722903, -1.0501579908443834,
								0.552104696331739, 1.3235049876698253, 1.653995689550814, -3.394333424870921, 
								-0.9610105412633638, 0.8093948807826427, 0.11769807478246722, -1.9790368431325849, 
								-0.7029974119950451, -1.4983149955801878, 0.6534878100392657, -3.891309775284998, 
								1.3222680744732622, 1.9183099725018982, -0.3929270788239533, -0.14594592899955827, 
								-0.5186136376055442, 0.9001027718643764, -0.5974427601607452, 0.2796801114496016, 
								-0.09391272946230618, -1.045242211961467, 2.593557991444645, -0.090813930607678, 
								1.0250362501432961, -0.17122832740598104, -1.198700482473792, -0.9809567820478314, 
								-1.8191247868912985, 0.8889126133184919, -0.6741655645323877, -1.693234113265235, 
								-0.49958698863835777, 0.7998022612825357, -3.0317911174199046, 1.1304942391301331, 
								-0.06961436942789863, -0.4310418205571006, 2.653037506118875, -1.1318004061135787, 
								-0.4114424920167265, -1.69383312062336, -0.01410251893470846, -0.5452812818663919, 
								1.0471150011773218, 1.0592070241592126, -0.7589476721645105, 0.18504815081645362, 
								-1.2782371758903799, -3.198580873383889, -1.1223973154546125, -2.139907402106115, 
								-1.454910257786402, 1.2396007359617567, -0.17099628216261573, -0.7149533883090601, 
								-0.5706033818418597, -0.7246795851955518, -0.03912070341748527, -0.7949524332954405,
								0.37773239941487924, 2.574906449233901, 0.827793259912873, -2.0309145745250032, 
								0.5575868762090828, -0.26688055267696165, 0.5998921304871374, 3.5363943360164614, 
								-0.9358474271885415, -2.7559661434208658, -2.1764407740474345, -4.258955343429767, 
								1.4096700130243005, 0.570221115797754, -1.0127710090178867, -1.5087202877102261, 
								1.0277432267859326, -0.00301073442000724, 1.3700734243788995, 1.1174775813935005, 
								2.7490397650854677, -0.14129260606422017, 1.1858708005085623, -2.66350892883937, 
								0.7970079009999912, 0.8398651635438542, 0.8459919232508679, 2.1023346075015326, 
								-0.20225230631943708, 0.2231229004566983, 0.3671199073910663, 1.2118463670962214, 
								-2.2384960504883185, 0.42230189718143096, -0.12773589214620915, -3.3737229340634274, 
								0.6712104822125662, -0.5623165862496264, 0.8247037699175153, -1.4888431086655922, 
								1.8893791827770754, 0.18355794578997292, 2.2260881572471645, -0.9858956590938109, 
								1.8746705867711673, -0.4492582324508535, -1.1353284228836888, 0.8165900279107424, 
								-0.6529368831573114, 0.3903348624914426, 2.2600258364623755, 2.0470077735087027, 
								-3.6950221827833567, -1.1783627598924027, -0.9744685016659325, -1.5994141130282433, 
								-1.6067067987309736, -0.9234047460975171, 0.8904211381603742, 0.36649576704799497, 
								-0.8908543097159564, 0.006751968596901002, -0.8280730163161361, 0.053153381596026894, 
								-1.6358883865084257, 1.775185154367792, -1.1352265177125673, 1.0979446235503794, 
								1.593879466651264, -0.598538569715115, 1.8458876019007844, -1.5119354046080975, 
								-2.020085731316158, -1.5820950070300517, 1.5921658487034873, -1.6915162486433535, 
								-0.8776761535967024, 0.2905131614488509, 1.3030356701078532, 2.0262747014717513, 
								0.7234069533735158, -1.3389030758657894, -0.02283449566532983, -0.839808852588597, 
								-0.6955680374222301, -0.41666363436016796, -0.7269541261448587, -0.5823206407143255, 
								-0.6117844621961339, -1.773583091824382, -1.8227078809643344, 1.4207520656307255, 
								1.0121021710607492, -1.0511076716718846, -0.8328616820090229, 0.5525231581006231, 
								0.22424468315478982, -0.12233810634368755, 0.8193805099463256, -0.5505150246836703, 
								0.35586119849503695, -0.37773817223006856, 0.6700310825563999, -0.25171926499957525, 
								-2.4650196915577953, 2.800026738632388, 0.054399205044297454, 1.344031410003991, 
								2.3992525017703126, -0.6144691768621632, -0.7980321405085672, 0.09294352462041401, 
								-1.5271191398217039, -1.8650809219275217, 1.4429701838881734, -0.422428624312937, 
								0.31407319277442136, 2.0669638157479717, -2.015618039025254, -0.13445688962889887, 
								1.3482233937754105, -0.3194690793857528, 0.21732613915348684, 2.442306142045435, 
								-1.5955875072461467, 1.1493621229640454, -0.5985399646708643, -2.2313444192462706, 
								0.001580597933558622, -1.5097576799709784, 0.055679886350673564, 1.349812668359207, 
								-0.3626116666139402, 1.2070011589841751, 1.3061669665122724, -0.5815901122746863, 
								0.21351275509397952, -0.675725729811047, -0.31106945219796817, -0.981390859401946, 
								-3.5857003244344434, -2.3794932715309702, 1.1410095362772694, -0.5061830578454546, 
								0.6247891279211158, -0.13929845419223807, 1.0095007205309476, -1.4498203255861808, 
								-1.2380237058087622, -1.3975512673926584, 3.6192686953475905, -0.8521422680139781, 
								-1.256510390904777, 0.6035364577595314, -1.897646611462231, 0.28906289437294486, 
								1.1044495340658074, 0.553553532127761, -4.164151939414819, -2.624226855698873, 
								1.1975605880750069, -0.7152829304252047, -2.708421507684577, -0.5510383273432158, 
								4.250844787135506, 1.1224862502644744, -1.13279081156131, 2.6099108976865266, 
								-0.6097371203551685, 0.8619095985000531, 1.2300165570396098, 0.6098789457204135, 
								0.3769957385970449, 1.1666204469249135, 0.11527274227211744, -1.7239920819564252, 
								-1.7926566184263342, 0.9490372958673449, 2.1190284957412273, -0.23938025221728576, 
								1.3428534781977988, 1.0717431372940582, -0.5310526699213246, -0.22330286448334508, 
								-1.6855723644039675, -1.180771513338842, -0.7487353752023637, -1.202993762676471, 
								1.1335731161210185, 0.23002972954548495, -1.3331189323123633, -0.4011973019224974, 
								0.8943790751048741, 0.515048962766933, -0.8829598686263169, -1.4006252379868074, 
								-2.5880406352361787, 0.1471094248703182, 1.1767842630861696, -1.5367586202119956, 
								0.7956207548603502, -1.5675189950039516, 0.47140515042798614, -0.5681143837173736, 
								-1.821821908478681, -0.6823417090901566, -0.8541361406907755, -1.1302245799877726, 
								0.7294090393147691, 2.876357134201796, -0.9484014901395184, 1.6957089184229452, 
								-0.12464592047460789, 1.7782771657641674, 1.6566005571247686, 1.4677115074761227, 
								0.31642335100266894, 0.3820080814495714, -0.041818008824671515, -1.257385962691819, 
								-0.06050303361865433, -0.3476519412806997, 1.8370485942451165, -0.5273778040565024, 
								1.7790469476991282, 1.2823436645930049, 1.1631736352488413, -0.4711068911658804, 
								1.2298476801025, -1.8353007215004595, -1.5246583553372721, -0.4566826998451468, 
								-0.2899587130110082, -1.9903729679687263, 0.43978662480277186, 0.5400074036265661, 
								-1.2013887842861999, -1.1455974060279976, -1.7444267261134758, -0.3570628568719624, 
								-1.7500012810313674, -0.6034232949152541, 0.05185351312676353, -0.15357622082078579, 
								-0.2273063078991938, -0.14252169083427899, 3.419085796745461, -1.1450640971157493, 
								-1.1077192722737819, 1.730977508585181, 0.8721324858910471, 1.4299834646299279, 
								-0.4102445411576849, -1.9159611680137156, 0.6984616025006694, 0.6543726189067457, 
								0.8570203204183897, 0.8385540877334238, 0.43228597244020617, 0.3537275969148874, 
								0.8747575625897029, 1.0934574284863119, -0.4776819142706088, -1.7818897586590603, 
								0.43877801571651104, 0.4735271288798478, -0.8167154345520726, -0.8843514485870845, 
								-3.302963777109195, -0.3980319861346514, 0.3444399697225685, -2.160854080065189, 
								-0.014248780175530882, -0.9001254868735152, 0.9815788953527007, -0.5503402338278136, 
								-1.4101936168539233, 0.6855579500367065, -0.08379671776065473, 0.7462712493176826, 
								0.0037434713519080435, 0.05993187683615106, -0.8102951976361906, -0.29041934207414716, 
								0.9503084455432103, 1.055376087234944, -0.9733375381242768, 0.02275160195636818, 
								0.18917222536855088, 1.440645790191831, 1.231252854850646, -1.650116469369853, 
								-0.23498532953074427, 1.6558322979686477, -1.9964456397067059, -0.6328433967931103, 
								-2.307212521578366, 0.7736956438080307, 0.16912976179267453, -0.2245174570103416, 
								0.43685031065127483, -0.19954969704903747, 2.7590491907795727, -1.5315850556053032, 
								0.14644923473924756, 1.5062562732087874, -0.42620193305663157, 0.21023973742195368, 
								-2.547725197698647, -0.6003118029142457, 0.8724193722340179, 1.0747321312800968, 
								-0.5590275430076012
		};		
		return MatrixOperations.convArrayToVec(whiteNoise);
	}
	
	
	public static double [][] sawtoothWave(int length, int period, double amplitude, double yShift) {
		double [][] sawtoothWave = new double [length][1];
		double stepFrac = amplitude/((double) period);
		int counter = 0;
		for(int i=0; i<length; i++) {					
			sawtoothWave[i][0] = (counter*stepFrac)+yShift;
			if(counter == period) {
				counter = 0;
			}
			counter++;
		}
		return sawtoothWave;
	}

	
	public static double [][] tangensWave(int length, int period, double amplitude, double yShift) {
		double min = -0.45*Math.PI;
		double max = 0.45*Math.PI;
		double s   = Math.tan(max);
		double diff = max-min;
		double step = diff/((double) period);
		double [][] tangensWave = new double [length][1];
		int counter = 0;
		for(int i=0; i<length; i++) {
			double x = min+counter*step;
			tangensWave[i][0] = Math.tan(x)/s*amplitude+yShift;		
			if(counter == period) {
				counter = -1;
			}
			counter++;
		}
		return tangensWave;
	}
	
	
	public static double [][] sinusWave(int length, int period, double amplitude, double yShift) {
		double min = 0.0;
		double max = 2.0*Math.PI;
		double step = max/((double) period);
		double [][] sinusWave = new double [length][1];
		int counter = 0;
		for(int i=0; i<length; i++) {
			double x = min+counter*step;
			sinusWave[i][0] = Math.sin(x)*amplitude+yShift;		
			if(counter == period) {
				counter = 0;
			}
			counter++;
		}
		return sinusWave;
	}
	
	
	public static void plot_ica(ArrayList<double [][]> factors, ArrayList<double [][]> reconstruction) {
		
		ArrayList<String []> titleList = new ArrayList<String []>();
		String [] titles = {"ICA", "ICA"};
		String [] subTitles = {"Extracted Factors", "Reconstructed Observations"};
		
		titleList.add(titles);
		titleList.add(subTitles);
		
		plotSignals(factors, reconstruction, titleList);	
	}
	
	
	public static void plot_pca(ArrayList<double [][]> factors, ArrayList<double [][]> reconstruction) {
		
		ArrayList<String []> titleList = new ArrayList<String []>();
		String [] titles = {"PCA", "PCA"};
		String [] subTitles = {"Extracted Factors", "Reconstructed Observations"};
		
		titleList.add(titles);
		titleList.add(subTitles);
		
		plotSignals(factors, reconstruction, titleList);
	}
	
	
	@SuppressWarnings("static-access")
	public static void plotSignals(ArrayList<double [][]> orgSignals, ArrayList<double [][]> obsSignals, ArrayList<String []> titleList) {
		
		GenGraphics graph = new GenGraphics();
		
		int nSignals = orgSignals.size();
		int nSamples = 2*nSignals;
      		
	 	String [] yLabels   = new String [nSamples];
	 	String [] xLabels   = new String [nSamples];
	 	
	 	graph.setNumberOfPlotColums(2);
	 	graph.setNumberOfPlotRows(nSignals);
	 	
	 	graph.setGraphWidth(700);
	 	graph.setGraphHeight(500);
	 	List<Color> defColors = graph.getDefaultColorsAsList();
	 	
	 	int n = orgSignals.get(0).length;
	 	double [][] x = graph.get_default_x_axis_labels(n);
	 	int counter = 0;
	 	int idx = 0;
	 	for(int i=0; i<nSamples; i++) {
	 		
	 		xLabels[i]   = "";
	 		yLabels[i]   = "";
		 	 
	 		double [][] y = new double [1][1];
	 		
	 		if(counter == 0) {
	 			y = orgSignals.get(idx);
	 		}
	 		if(counter == 1) {
	 			y = obsSignals.get(idx);
	 		}
	 				 		
	 		graph.plotLines(x, y, true, defColors.get(idx));
	 		
	 		if(counter == 1) {
	 			idx++;
	 			counter = -1;
	 		}
	 		
	 		counter++;
	 	}
	 	
	 	String [] titles = new String [2];
	 	String [] subTitles = new String [2];
	 	
	 	if(titleList == null) {
	 		titles[0] = "Source Signals";
		 	titles[1] = "Observed Signals";
	 	}else {
	 		if(titleList.get(0) == null) {
		 		titles[0] = "Source Signals";
			 	titles[1] = "Observed Signals";
	 		}else {
		 		titles[0] = titleList.get(0)[0];
			 	titles[1] = titleList.get(0)[1];
	 		}
	 	}
	 	if(titleList == null) {
	 		subTitles[0] = "";
	 		subTitles[1] = "";	
	 	}else {
	 		if(titleList.get(1) == null) {
		 		subTitles[0] = "";
		 		subTitles[1] = "";	
	 		}else {
		 		subTitles[0] = titleList.get(1)[0];
		 		subTitles[1] = titleList.get(1)[1];	
	 		}
	 	}
	 	
	 	
	 	//graph.setLineColor(lineColor);	 	
	 	graph.setBackroundColor(Color.WHITE);
	 	graph.setTitle(titles, null, "11");
	 	graph.setSubTitle1(subTitles, null, "10");
	 	graph.setYLabel(yLabels, null, "8");
	 	graph.setXLabel(xLabels, null, "8");
	 	graph.setNumberOfDigits4XAxis(0);   
	 	graph.setNumberOfDigits4YAxis(1);
	 	graph.setFontOfXAxisUnits("plain", 8);
	 	graph.setFontOfYAxisUnits("plain", 8);
	 	graph.set_point_width(2);
	 	graph.set_line_widht(1);
	 	
	 	graph.plot();
	 		 	
	}
	
	
	@SuppressWarnings("static-access")
	public static void plotDistribution4ICA(double [][] x, double [][] distValues, String distribution) {
		
		GenGraphics graph = new GenGraphics();		
	 	graph.setNumberOfPlotColums(1);
	 	graph.setNumberOfPlotRows(1); 	
	 	graph.setGraphWidth(700);
	 	graph.setGraphHeight(500);
	 	
	 	graph.plotLines(x, distValues, true, Color.BLUE);
	 		
	 	String titles [] = new String [1];
	 	String yLabels [] = new String [1];
	 	String xLabels [] = new String [1];
	 	titles[0] = distribution;
	 	yLabels[0] = distribution;
	 	xLabels[0] = "";
	 			
	 	//graph.setLineColor(lineColor);	 	
	 	graph.setBackroundColor(Color.WHITE);
	 	graph.setTitle(titles, null, "11");
	 	graph.setYLabel(yLabels, null, "8");
	 	graph.setXLabel(xLabels, null, "8");
	 	graph.setNumberOfDigits4XAxis(0);   
	 	graph.setNumberOfDigits4YAxis(1);
	 	graph.setFontOfXAxisUnits("plain", 8);
	 	graph.setFontOfYAxisUnits("plain", 8);
	 	graph.set_point_width(2);
	 	graph.set_line_widht(1);
	 	
	 	graph.plot();
	 		 	
	}
	
	
	public static void main(String[] args) throws Exception {		
		test1();
		//test2();
	}
	
}

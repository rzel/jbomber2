
class Format142 {
	public static final String strformat(String s, int[] iarr) {
		Object oarr[] = new Object[iarr.length];
		for (int i = 0; i < iarr.length; ++i)
			oarr[i] = new Integer(iarr[i]);
		return String.format(s, oarr);
	}

	public static final String strformat(String s, double[] iarr) {
		Object oarr[] = new Object[iarr.length];
		for (int i = 0; i < iarr.length; ++i)
			oarr[i] = new Double(iarr[i]);
		return String.format(s, oarr);
	}

	public static final void strprintf(String s, int[] iarr) {
		s = strformat(s, iarr);
		System.out.print(s);
	}
}
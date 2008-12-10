import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.border.*;
import javax.swing.table.*;

public class IsoLineSelectPanel extends ColormapSelectPanel {
	private UnboundedDoubleSpinnerModel minIsoValueSpinnerModel = new UnboundedDoubleSpinnerModel(0.0);
	private UnboundedDoubleSpinnerModel maxIsoValueSpinnerModel = new UnboundedDoubleSpinnerModel(1.0);
	private SpinnerNumberModel isoLineCountSpinnerModel         = new SpinnerNumberModel(6,1,100,1);

	public IsoLineSelectPanel(int minColor, int maxColor, int colorCount, int colormap, JFrame frame) {
		super(minColor,maxColor,colorCount,colormap,frame);
		add(initIsoValuesSelectPanel());
	}

	public JPanel initIsoValuesSelectPanel() {
		JPanel isoValueSelectPanel = new JPanel();
		isoValueSelectPanel.setLayout(new GridLayout(3, 2));

		isoValueSelectPanel.add(new JLabel("Low: ", SwingConstants.RIGHT));
		isoValueSelectPanel.add(new JSpinner(minIsoValueSpinnerModel));
		minIsoValueSpinnerModel.setValue(new Double(0.0));
		isoValueSelectPanel.add(new JLabel("High: ", SwingConstants.RIGHT));
		isoValueSelectPanel.add(new JSpinner(maxIsoValueSpinnerModel));
		maxIsoValueSpinnerModel.setValue(new Double(1.0));


		/****************************
		 ** BEGIN Default gradient **
		 ****************************/
		DefaultTableModel tablemodelcolors = (DefaultTableModel)colortable.getModel();
		TableColumnModel cm = colortable.getColumnModel();
		ColorTableRenderer r = new ColorTableRenderer();
		TableColumn c = cm.getColumn(0);
		c.setCellRenderer(r);

		tablemodelcolors.setRowCount(0);
		Object[] o = new Object[1];
		o[0] = new Color(255, 0, 0);
		tablemodelcolors.addRow(o);
		o[0] = new Color(255, 255, 255);
		tablemodelcolors.addRow(o);
		colorCountSlider.setMinimum(colortable.getRowCount() - 1);

		//generate_custom_gradient_cache(); // handled by setInterpolateMode
		setInterpolateMode(ColormapSelectPanel.INTERPOLATE_HSV);
		/**************************
		 ** END Default gradient **
		 **************************/


		isoValueSelectPanel.add(new JLabel("Choose n:", SwingConstants.RIGHT));
		isoValueSelectPanel.add(new JSpinner(isoLineCountSpinnerModel));

		isoValueSelectPanel.setBorder(new TitledBorder("Iso values:"));
		isoValueSelectPanel.setAlignmentX(Component.LEFT_ALIGNMENT);
		return isoValueSelectPanel;
	}

	public double getMinIsoValue() {
		return ((Double)minIsoValueSpinnerModel.getValue()).doubleValue();
	}

	public double getMaxIsoValue() {
		return ((Double)maxIsoValueSpinnerModel.getValue()).doubleValue();
	}

	public int getIsoLineCount() {
		return ((Integer)isoLineCountSpinnerModel.getValue()).intValue();
	}
}

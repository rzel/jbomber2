import java.io.*;
import java.net.*;

import javax.swing.*;

import java.awt.*;
import java.awt.event.*;

/**
 * A <code>GameComponent</code> is an object that is used to
 * draw the game. A <code>KeyListener</code> is specified here, to
 * handle player input.
 */
public class GameComponent extends JComponent implements Runnable, KeyListener
{
    // Something
    public static final long serialVersionUID = 1;

    // The Tread to play the game
    private Thread game_thread;

    // The DelegatorClient that is used to send data to server
    private DelegatorClient delegator_client;

    /**
     * Creates a new <code>GameComponent</code> object.
     */
    public GameComponent(DelegatorClient delegator_client)
    {
        this.delegator_client = delegator_client;
        setFocusable(true);
        addKeyListener(this);
        start();
    }

    /**
     * Start the thread to play this game.
     */
    public void start()
    {
        game_thread = new Thread(this);
        game_thread.start();
    }

    /**
     * Draw the game
     */
    public void run()
    {
        for (;;)
        {
            repaint();

            try
            {
                game_thread.sleep(20);
            }
            catch(InterruptedException ie) {}
        }
    }

    /**
     * Handles players' keypresses
     */
    public void keyPressed(KeyEvent e)
    {
        // Only handle data when it is possible
        if (!delegator_client.canHandleData()) return;

        switch (e.getKeyCode())
        {
            case Config.MOVE_RIGHT:
                delegator_client.handleAction(Config.ACT_MOVE_RIGHT);
                break;
            case Config.MOVE_LEFT:
                delegator_client.handleAction(Config.ACT_MOVE_LEFT);
                break;
            case Config.MOVE_DOWN:
                delegator_client.handleAction(Config.ACT_MOVE_DOWN);
                break;
            case Config.MOVE_UP:
                delegator_client.handleAction(Config.ACT_MOVE_UP);
                break;
        }
    }

    /**
     * Draws the scheme used in this game.
     *
     * @param g The <code>Graphics</code> to draw on.
     */
    private void drawScheme(Graphics g)
    {
        // TODO good drawing
        char[][] scheme = delegator_client.getScheme();

        if (scheme == null) return;

        for(int y = 0; y < scheme.length; ++y)
        {
            for(int x = 0; x < scheme[y].length; ++x)
            {
                switch (scheme[y][x])
                {
                    case Scheme.BRICK_CHAR:
                        g.setColor(Color.GRAY);
                        break;
                    case Scheme.SOLID_CHAR:
                        g.setColor(Color.DARK_GRAY);
                        break;
                    default:
                        g.setColor(Color.GREEN);
                }

                g.fillRect(y * Config.GRID_WIDTH,
                           x * Config.GRID_HEIGHT,
                           Config.GRID_WIDTH,
                           Config.GRID_HEIGHT);

            }
        }
    }

    /**
     * Draws a player.
     *
     * @param g The <code>Graphics</code> to draw on.
     * @param p The <code>ReadablePlayer</code> to draw.
     */
    private void drawPlayer(Graphics g, ReadablePlayer p)
    {
        // TODO good drawing
        g.setColor(Color.BLUE);
        g.fillOval(p.getXPos(),
                   p.getYPos(),
                   Config.GRID_WIDTH,
                   Config.GRID_HEIGHT);

        g.drawString(p.getNickname(), p.getXPos(), p.getYPos());
    }

    /**
     * Paint the game
     *
     * @param g The graphics to paint the game to.
     */
    public void paint(Graphics g)
    {
        g.setColor(Color.RED);
        g.fillRect(0,0,(int)getSize().getWidth(), (int)getSize().getHeight());

        if (delegator_client.getScheme() != null)
        {
            drawScheme(g);
        }
        else
        {
            g.drawString("Waiting...", 100, 100);
        }

        for (ReadablePlayer p : delegator_client.getPlayers())
        {
            drawPlayer(g, p);
        }
    }

    /**
     * Not implemented
     */
    public void keyTyped(KeyEvent e) {}


    /**
     * Not implemented
     */
    public void keyReleased(KeyEvent e) {}

    //TESTING
    public static void main(String[] args)
    {
        DelegatorClient client;

        try
        {
            String host;

            if (args.length > 0)
            {
                host = args[0];
            }
            else
            {
                host = Config.SERVER_ADDR;
            }

            InetAddress addr = InetAddress.getByName(host);
            client = new DelegatorClient(addr,
                                         Config.SERVER_PORT,
                                         System.out,
                                         "Test");

            JFrame f = new JFrame("BOMBERMAN");
            f.setSize(500,500);
            f.setLocationRelativeTo(null);
            f.setLayout(null);

            GameComponent gc = new GameComponent(client);
            gc.setBounds(0,0,Config.GAME_WIDTH, Config.GAME_HEIGHT);
            f.add(gc);
            f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            f.setVisible(true);
        }
        catch (IOException ioe)
        {
            // If the server could not be started, give an error message.
            System.out.println(ioe.getMessage() + "\n");
            System.exit(-1);
        }
    }
}